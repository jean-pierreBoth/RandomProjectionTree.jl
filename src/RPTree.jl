#  File for random projection tree implementation
# we instantiate with Tree and node with T=Vector{Float64}

# interesting distance are
#                        Euclidean     L2
#                        Cityblock     L1
#                        Jaccard       fuzzy nucleotide : 1 - sum min (a,b)/ sum max(a,b) 
#                        CosineDist
# they are instantiated by calling Euclidean(),  ...
#
#  TODOs:
#     diameter estimation with kcover algorithm?
#     optimisation of threshold used in splitting
#


rptreeDebugLevel = 1
debuglock=Threads.RecursiveSpinLock()


include("Tree.jl")

using Distances
using Clustering
using Match
using Distributed
using Statistics
using LinearAlgebra
using Random


const splitDiam = Int64(1)
const splitProj = Int64(2)
const splitNull = Int64(3)   # when a node is not splitted (leaf)

export
RPTree,
RPTreeArg,
RPTevent,
RPTProjParams,
RPTreeEvent,
randomProjection,
analyzeSplittingInfo


"""
KeyVector

base data of node is an indexed vector
So it is possible when in a node to know the rank
in original data in tree the node consists in 
keys(KeyVector) consists in rank of data in the  node in the original vector of vectors 
dispatched in the tree  
"""
const KeyVector = Dict{Int64,Vector{Float64}}


"""

# RPTDiamSplit

This struct stores data necessary to describe a split by diameter

    FIELDS
    -----

* pivot : the vector aroud which we construct a ball to split
* radius : distance to pivot above which we go to right, under which we go to left

    """
mutable struct RPTDiamSplit
    pivot::Vector{Float64}
    radius::Float64
    #
    function RPTDiamSplit(center::Vector{Float64})
        new(center, -1.)
    end
end


"""
# RPTProjParams

This struct stores the dalong which data are projected.

 FIELDS
 -----

    * projray : ray on which do project data
    * cosineVector : a vector containing for each data the cosine with projray.

     Data are split according to position of cosine above or below median of cosine)
"""
mutable struct RPTProjParams
    projray::Vector{Float64}
    cosineVector::Vector{Float64}
    #
    function RPTProjParams(ray::Vector{Float64}, cosVec::Vector{Float64})
        new(ray,cosVec)
    end
end



"""

# RPTreeEvent

Event on Node of Tree
We store the event as it enables propagating a new data across the tree if we
want to do any learning

 FIELDS
 -----

* split : the split event coded as constants `splitDiam` or `splitProj` or `splitNull`
* diameters : a vector of size 2 containing mean diameter and then max diameter
* projEvent  field of struct Union{RPTProjParams, Nothing}
* diamEvent  field of struct Union{RPTDiamSplit, Nothing}

  METHODS
  ------

    1. `function RPTreeEvent(s::Int64, diameters::Vector{Float64}, pivot::Union{RPTDiamSplit,Nothing})`
    constructor for split on diameter criteria

    2. `function RPTreeEvent(s::Int64, diameters::Vector{Float64}, projData::Union{RPTProjParams, Nothing})`
        constructor for split on projection against a random Vector
        as diameters are always estimated (the struct of split depends upon computed diameters),
        they always are in the constructor

    3. `function RPTreeEvent(diameters::Vector{Float64})`
         a constructor used to just store diameter estimation at leaves. Used with event splitNull
"""
mutable struct RPTreeEvent
    split::Int64
    # first mean diam then max diam
    diameters::Vector{Float64}
    projEvent::Union{RPTProjParams, Nothing}
    diamEvent::Union{RPTDiamSplit, Nothing}
    # constructor for split by diameter
    function RPTreeEvent(s::Int64, diameters::Vector{Float64}, diamEvent::Union{RPTDiamSplit, Nothing})
        # could check that s is splitDiam or splitProj
        new(s, diameters, nothing, diamEvent)
    end
    #
    function RPTreeEvent(s::Int64, diameters::Vector{Float64}, projData::Union{RPTProjParams, Nothing})
        new(s, diameters, projData)
    end
    #
    function RPTreeEvent(diameters::Vector{Float64})
        new(splitNull, diameters, nothing, nothing)
    end
end

const  RPTNode = TreeNode{KeyVector, RPTreeEvent}


"""

# RPTreeArg
 collects the argument for tree growing

 FIELDS
 -----
 *  D : the metric to use to compute all distances
 *  depth : The maximum depth of the tree
 *  threshold : the ratio between maxDiameter and mean distance between 2 objects above which
    we split by diameter. If mean / maxDiameter > threshold we split by diameter. So
    with threshold = 1. we always use diameter rule and the higher threshold is more 
    we use projection

 CONSTRUCTORS
 -----------
 
  
"""
mutable struct RPTreeArg
    D::Distances.SemiMetric
    # the depth of tree
    depth::Int64
    # the ratio between maxDiameter and mean distance between 2 objects above which we split by diameter
    threshold::Float64
    #
    RPTreeArg(D::Distances.SemiMetric, depthArg::Int64, thresholdArg::Float64) = new(D, depthArg, thresholdArg)
end


"""
# RPTree
The struct RPTree stores the list of event that occurred during tree construction
The parameters used to grow the tree

 FIELDS
------    
* treedata : the Tree 
* argument : the parameters describing tree growing
* eventDict : the event that occured in each node.
* treelock : a lock to manipulate eventDict 
    
"""
mutable struct RPTree
    treedata::Tree{KeyVector, RPTreeEvent}
    argument::RPTreeArg
    eventDict::Dict{TreeNode{KeyVector,RPTreeEvent}, RPTreeEvent}
    #
    treelock::Threads.RecursiveSpinLock
    #
    function RPTree(argumentArg::RPTreeArg, vdata::Array{Vector{Float64},1})
        eventDict = Dict{TreeNode{KeyVector,RPTreeEvent}, RPTreeEvent}()
        # transform data into indexed data
        keyData = KeyVector()
        for i in 1:length(vdata)
            keyData[i] = vdata[i]
        end
        new(Tree{KeyVector, RPTreeEvent}(keyData), argumentArg, eventDict, Threads.RecursiveSpinLock())
    end
end



function myjaccard(a::Array{Float64,1}, b::Array{Float64,1})
    num = 0.
    den = 0.
    for I in 1:length(a)
        @inbounds ai = a[I]
        @inbounds bi = b[I]
        abs_m = abs(ai-bi)
        abs_p = abs(ai+bi)
        num += abs_p - abs_m
        den += abs_p + abs_m   
    end
    1. - num/den
end


function jaccard_C(a::Array{Float64,1}, b::Array{Float64,1})
    return ccall((:jaccard,"mylibjulia.so"),Float64,(Int64,Ptr{Float64},Ptr{Float64}),length(a),a,b)
end

#
#  estimation of node center and mean distance to center
#  This loop should be parallelized
#





#
#  
#



"""
diameterEstimation 

compute some  estimates of distances between pairs of element in node

returns
* mean distance between pairs,
* max distance of pair
* key of data element with min sum of distances to others.
"""
function diameterEstimation(D::Distances.SemiMetric, node::TreeNode{KeyVector} , fraction::Float64)
    #
    # diameter max and mean estimation
    # We sample elements and compute distances between all couples
    #                                          
    nb=length(node.data)
    if nb <= 1
        return 0.,0.,0
    end
    allKeys = collect(keys(node.data))
    keptI = filter(x -> (rand() < fraction) , allKeys)
    # could not use map!(node -> evaluate(D, node.data, center) , distances, collection)
    # as we do not use all indices. ??
    nbkept=length(keptI)
    costs=zeros(nbkept,nbkept)   # to store sum of distance of each kept index to all others.
    maxByI=zeros(nbkept)
    maxDiameter=0.
    meanDiameter=0.
    nbcouple = 0
    for i in 2:nbkept
        distancesToI = Vector{Float64}(undef,i-1)    # in fact this is cost free versus one maximal allocation
                                               # and maximum taken on a subrange
        map!(j->evaluate(D, node.data[keptI[i]], node.data[keptI[j]]), distancesToI, 1:i-1)
        costs[1:i-1,i]=distancesToI
    end
    Fcosts=Symmetric(costs,:U)
    maxDiameter = maximum(Fcosts)
    meanDiameter = sum(Fcosts)
    # as we have 0 on diagonal, we divide by full matrix minus diagonal
    nbcouple = (nbkept-1) * nbkept
    meanDiameter = meanDiameter/nbcouple
    costmin,imin = findmin(costs)
    #
    return meanDiameter, maxDiameter , keptI[imin] 
end



# the following is faster than the preceding as soon as we have 1000 nodes after applying fraction
# for vectors of lenght 10000. So we need 10**5 distances to use parallelism
#

# in the following we should use remotecall_fetch with a

function diameterEstimationByRange(D::Distances.SemiMetric, node::TreeNode{KeyVector,RPTreeEvent},
                                   keptI::Array{Int64,1}, range::UnitRange{Int64})
    maxByTask = 0.
    sumByTask = 0.
    nbkept = length(keptI)
    costs=zeros(nbkept,nbkept)   # to store sum of distance of each kept index to all others.
    ikeys = keys(node.data)
    for i in range
        distancesToI = Vector{Float64}(undef, i-1)
        map!(j->evaluate(D, node.data[ikeys[keptI[i]]], node.data[ikeys[keptI[j]]]), distancesToI, 1:i-1)
        costs[1:i-1,i]=distancesToI
    end
    Fcosts=Symmetric(costs,:U)
    maxByTask = maximum(Fcosts)
    # sum and convert the 1*nbkept Array{Float64,2} to a Vector{Float64} !!
    colCost=convert(Vector{Float64},sum(Fcosts,1)[1,:])
    sumByTask = sum(colCost)
    #
    return sumByTask, maxByTask, colCost
end


# just to test remote call , not really efficient

function diameterEstimation2tasks(D::Distances.SemiMetric, node::TreeNode{KeyVector,RPTreeEvent} , fraction::Float64)
    #
    # diameter max and mean estimation
    # We sample elements and compute distances between all couples
    #                                          
    nb=length(node.data)
    keptI = filter(x -> (rand() < fraction) , 1:nb)
    # could not use map!(node -> evaluate(D, node.data, center) , distances, collection)
    # as we do not use all indices. ??
    nbkept=length(keptI)
    maxDiameter=0.
    meanDiameter=0.
    nbcouple = 0
    res = Vector{ Tuple{Float64,Float64, Vector{Float64}} }(undef, 2)
    m=round(Int64, nbkept/sqrt(2))
    @sync begin
        @async  res[1] = remotecall_fetch(2, diameterEstimationByRange, D, node, keptI, 2:m)
        @async  res[2] = remotecall_fetch(3, diameterEstimationByRange, D, node, keptI, m+1:nbkept)        
    end
    meanDiameter = res[1][1] + res[2][1]        # sum of means
    maxDiameter = max(res[1][2], res[2][2])     # max of max
    nbcouple = (nbkept-1) * nbkept
    meanDiameter = meanDiameter/nbcouple
    # compute costs, get full symmetric matrix
    costs=res[1][3] .+ res[2][3]
    costmin,imin = findmin(costs)    
    #
    return meanDiameter, maxDiameter, imin
end

#
# The following is not efficient
#

function diameterEstimation2Threads(D::Distances.SemiMetric, node::TreeNode{KeyVector,RPTreeEvent} , fraction::Float64)
    #
    # diameter max and mean estimation
    # We sample elements and compute distances between all couples
    #                                          
    nb=length(node.data)
    keptI = filter(x -> (rand() < fraction) , 1:nb)
    # could not use map!(node -> evaluate(D, node.data, center) , distances, collection)
    # as we do not use all indices. ??
    nbkept=length(keptI)
    maxDiameter=0.
    meanDiameter=0.
    nbcouple = 0
    res = Vector{ Tuple{Float64,Float64, Vector{Float64}} }(undef, 2)
    # compute m so that the 2 threads have the same amount of work!!!
    m=round(Int64, nbkept/sqrt(2.))
    myRanges=Array{UnitRange{Int64}}(2)
    myRanges[1] = 2:m
    myRanges[2] = m+1:nbkept
    Threads.@threads  for i in 1:2  
        res[1] = diameterEstimationByRange(D, node, keptI, 2:m)
        res[2] = diameterEstimationByRange(D, node, keptI, m+1:nbkept)        
    end
    meanDiameter = res[1][1] + res[2][1]        # sum of means
    maxDiameter = max(res[1][2], res[2][2])     # max of max
    nbcouple = (nbkept-1) * nbkept
    meanDiameter = meanDiameter/nbcouple
    # compute costs
    costs=res[1][3] .+ res[2][3]
    costmin,imin = findmin(costs)    
    #
    return meanDiameter, maxDiameter, imin
end



function splitNodeByDiameter(D::SemiMetric, node::TreeNode{KeyVector, RPTreeEvent}, pivot::Int64)
    #
    #
    # compute mean
    #
    nb=length(node.data)
    center=node.data[pivot]
#    center=node.data[1]
#    for i in 2:nb
#        center = center .+ node.data[i] 
#    end
#    center = center ./ nb
    #
    # compute median to distances to center
    #
    distances = zeros(Float64,nb)
    vectors = collect(values(node.data))
    map!(vec -> evaluate(D, vec, center) , distances, vectors)
    medianDist=median(distances)
    #
    leftData = KeyVector()
    rightData = KeyVector()
    # it is documented that keys and values are in the same order !!
    # So the i-th distance correspond to the i-th values and to the i-th key
    vindex = collect(keys(node.data))
    for i in 1:length(vindex)
        if distances[i] <= medianDist
            # go to left
            leftData[vindex[i]] = vectors[i]
        else
            # go to right
            rightData[vindex[i]] = vectors[i]
        end
    end
    #
    # We can create children
    #
    leftNode = TreeNode{KeyVector, RPTreeEvent}(node, leftData)
    rightNode = TreeNode{KeyVector, RPTreeEvent}(node, rightData)
    #
    node.data = KeyVector()   # free .... all data have been pushed down
    #
    return leftNode,rightNode,medianDist
end # end of splitByDiameter



# generate random direction too split data in node

function generateRandomDirection2(node::TreeNode{KeyVector})
    item, state = iterate(node.data)    
    dim = length(item[2])
    #
    direction = zeros(dim)
    randomNumbers = Vector{Float64}(undef, dim+1)
    #
    Random.rand!(randomNumbers)
    i = 1
    while i < dim
        s = sqrt(-2 * log(randomNumbers[i]))
        t = 2*randomNumbers[i+1]*pi
        g1 = s * cos(t)
        g2 = s * sin(t)
        direction[i] = g1
        if i+1 <= dim
            direction[i+1] = g2
            i = i+2
        else
            break
        end
    end    # end while
    return direction
end




function generateRandomDirection(node::TreeNode{KeyVector, RPTreeEvent})
    item, state = iterate(node.data)    
    dim = length(item[2])
    #
    direction = zeros(dim)
    randomNumbers = Vector{Float64}(undef, dim+1)
    #
    Random.rand!(randomNumbers)
    halfsize = round(Int64,dim/2)
    @simd for i in 1:halfsize
        j=2*i-1
        @inbounds s = sqrt(-2 * log(randomNumbers[j]))
        @inbounds t = 2*randomNumbers[j+1]*pi
        g1 = s * cos(t)
        g2 = s * sin(t)
        @inbounds direction[j] = g1
        @inbounds direction[j+1] = g2        
    end
    if dim%2 != 0
        s = sqrt(-2 * log(randomNumbers[dim]))
        t = 2*randomNumbers[dim+1]*pi
        g1 = s * cos(t)
        direction[dim] = g1
    end
#
    return direction
end



function splitNodeByProjection(node::TreeNode{KeyVector, RPTreeEvent})
    # generate random normal direction
    direction = generateRandomDirection(node)
    # compute scalar product of each element in node with direction. we have the dot function
    nb=length(node.data)
    scalarProduct = Vector{Float64}(undef, nb)
    # possibly iterate to avoid the collect ?
    map!(x-> dot(x,direction), scalarProduct, collect(values(node.data)))                            
    medVal=median(scalarProduct)
    # dispatch
    leftData = KeyVector()
    rightData = KeyVector()
    indexes = collect(keys(node.data))
    #
    projEvent = RPTProjParams(direction,scalarProduct) 
    #
    for i in 1:nb
        if scalarProduct[i] <= medVal
            # go to left
            leftData[indexes[i]] = node.data[indexes[i]]
        else
            # go to right
            rightData[indexes[i]] = node.data[indexes[i]]
        end
    end
    #
    # We can create children
    #
    leftNode = TreeNode{KeyVector, RPTreeEvent}(node, leftData)
    rightNode = TreeNode{KeyVector, RPTreeEvent}(node, rightData)
    #
    #
    node.data = KeyVector()   # free .... all data have been pushed down
    #
    return leftNode,rightNode,projEvent 
end

 
# splitNodeDiamAndProjection
# Here we must choose the way we do the splitting


function splitNodeDiamAndProjection(rptarg::RPTreeArg , node::TreeNode{KeyVector})
    #
    D = rptarg.D
    threshold = rptarg.threshold
    split = splitNull
    private = nothing
    #
    if length(node.data) == 0
        exit(1)
    end
    #  do an estimation of diameter
    sampleSize = 512.
    fraction = min(1., sampleSize/length(node.data))
    #
#    medDiam, maxDiam, pivot  =  diameterEstimationThreaded(D, node, fraction)
#     medDiam, maxDiam, pivot  =  diameterEstimation2Threads(D, node, fraction)
    medDiam, maxDiam, pivot  = diameterEstimation(D, node, fraction)
    #
    # to store node properties
    #
    diameters=Vector{Float64}(undef, 2)
    diameters[1] = medDiam
    diameters[2] = maxDiam
    # depending on node sphericity we split one way or the other
    if maxDiam/medDiam >= threshold
        # as splitNodeByDiameter clears its data we must get at pivot before the split
        diamEvent = RPTDiamSplit(node.data[pivot])
        leftNode,rightNode, radius = splitNodeByDiameter(D, node, pivot)
        diamEvent.radius = radius
        split = splitDiam
        private = RPTreeEvent(split , diameters, diamEvent)
    elseif medDiam > 0
        leftNode,rightNode, projparams = splitNodeByProjection(node)
        split = splitProj
        private = RPTreeEvent(split , diameters, projparams)
        # else nothing , means all element are equal
    end
    # 
    node.private = private   
    #
    return leftNode,rightNode
end # end of splitNodeDiamAndProjection




# once we got in serial mode we stick to it

function splitNodeSerial(rptarg::RPTreeArg, node::TreeNode)
    #
    depth = rptarg.depth
    if rptreeDebugLevel > 0
        lock(debuglock)
        @printf stdout "\n splitNodeSerial process id : %d thread %d depth %d " myid() Threads.threadid()  node.depth
        unlock(debuglock)
    end
    #
    if node.depth < depth && length(node.data) > 0
        leftNode,rightNode = splitNodeDiamAndProjection(rptarg ,node)
        # possibly a lock on tree , but we have a lock on node in which we insert
        insertChildInNode(node, leftNode)
        insertChildInNode(node, rightNode)
        #
        splitNodeSerial(rptarg, leftNode) 
        splitNodeSerial(rptarg, rightNode)
    end
    if rptreeDebugLevel > 0
        lock(debuglock)
        @printf stdout "\n splitNodeSerial process id : %d thread %d depth %d end " myid() Threads.threadid()  node.depth
        unlock(debuglock)
    end
    return node
end



# The following can replace splitNodeSerial , it just do one step by hand and call 2 threads for
# further splitting
#
function splitNodeThreaded(rptarg::RPTreeArg, node::TreeNode)
    depth = rptarg.depth
    if rptreeDebugLevel > 0
        lock(debuglock)
        @printf stdout "\n splitNodeTheaded process id : %d depth %d " myid()   node.depth
        unlock(debuglock)
    end
    #
    if node.depth < min(depth, 3) &&  Threads.nthreads() > 1
        leftNode,rightNode = splitNodeDiamAndProjection(rptarg ,node)
        # possibly a lock on tree , but we have a lock on node in which we insert
        insertChildInNode(node, leftNode)
        insertChildInNode(node, rightNode)
        # now a threaded loop
        nodesToSplit = Vector{TreeNode}(undef, 2)
        nodesToSplit[1] = leftNode
        nodesToSplit[2] = rightNode
        Threads.@threads for i in 1:2
            splitNodeThreaded(rptarg,nodesToSplit[i]) 
        end
    else
        splitNodeSerial(rptarg, node)
    end
    if rptreeDebugLevel > 0
        lock(debuglock)
        @printf stdout "\n splitNodeTheaded process id : %d depth %d end " myid()   node.depth
        unlock(debuglock)
    end
    return node
end



function splitNodeParallel(rptarg::RPTreeArg, node::TreeNode)
    #
    depth = rptarg.depth    
    @printf stdout "\n  splitNodeParallel process id : %d depth %d " myid()  node.depth
    #
    if node.depth < depth && length(node.data) > 0
        leftNode,rightNode = splitNodeDiamAndProjection(rptarg, node)
        #
        if nworkers() >= 2 && node.depth < 1
            @sync begin
                resLeftSubTree =  @spawn splitNodeThreaded(rptarg, leftNode) 
                resRightSubTree = @spawn splitNodeThreaded(rptarg, rightNode)
            end
            leftNode2 = fetch(resLeftSubTree)
            rightNode2 = fetch(resRightSubTree)
            insertChildInNode(node, leftNode2)
            insertChildInNode(node, rightNode2)
        else
            insertChildInNode(node, leftNode)
            insertChildInNode(node, rightNode)
            splitNodeSerial(rptree, leftNode)
            splitNodeSerial(rptree, rightNode)
        end
    end
    return node
end



#
# The driver method that grow a tree up to argument depth
#

"""

function randomProjection(rptree::RPTree)

The driver method that grow a tree up to argument depth

NOTA : The method cannot be run twice on the same rptree


"""
function randomProjection(rptree::RPTree)
    #
    node=rptree.treedata.root
    #
    if nworkers() <= 1 && Threads.nthreads() <= 1
        @printf stdout "\n going serial"
        splitNodeSerial(rptree.argument, node)
    elseif nworkers() <= 1 && Threads.nthreads() > 1
        splitNodeThreaded(rptree.argument, node)
    else
        node = splitNodeParallel(rptree.argument, node)
        rptree.treedata.root = node
    end
    # now we must iterate through leaves of tree, compute centers and send that to clustering
    @printf stdout "\n randomProjection, collecting leaves of tree \n"
    leafCenters=Array{Vector{Float64},1}()
    leaf = getFirstLeftLeaf(rptree.treedata)
    while leaf != nothing
#        @printf stdout "\n depth %d rank in Parent %d" leaf.depth leaf.rankInParent
        realdata = values(leaf.data)        
        center = sum(realdata)
        center = center ./ length(realdata)
        push!(leafCenters,center)
        leaf=getNextLeafRight(rptree.treedata, leaf)
    end
    @printf stdout "\n randomProjection, collected nb leaves = %d \n" length(leafCenters)
    #
    return leafCenters
end




function getFirstLeftNodeAtDepth(rptree::RPTree, depth::Int64)
    if depth > rptree.argument.depth
        return NullNode
    end
    #
    node = rptree.treedata.root
    while node.depth < depth
        if length(node.children) > 0
            node = node.children[1]
        else
            return NullNode
        end
    end
    if node.depth == depth
        return node
    end
    #
    return NullNode
end



# to get diameter of terminal leaves)
function getLeafStatistics(rptree::RPTree)
end

#  a traversal of tree after splitting to get all splitting info (diameters and split mode)


"""

    `function fillSplittingInfo(rptree::RPTree)`

fills `rptree.eventDict` with events describing how each splitted node was split
For leaves computes median and max diameter and associate an event with splitNull


"""

function fillSplittingInfo(rptree::RPTree)
    node = rptree.treedata.root
    sampleSize = 512
    #
    nbseen = 0
    nbproj = 0
    nbdiam = 0
    while node != nothing
        if rptreeDebugLevel > 1
            @printf stdout " \n dump de noeud : %d, depth %d"  nbseen  node.depth
            dumpPos(node)
        end
        nbseen += 1
        if length(node.children) > 0
            rptree.eventDict[node] = node.private
        else
            # we have a leaf
            fraction = min(1., sampleSize/length(node.data))
            res = diameterEstimation(rptree.argument.D, node, fraction)
            diameters=[res[1] , res[2]]
            # the following just store diameters but do not register an event so that analyzeSplittingInfo
            # will not be polluted
            event = RPTreeEvent(diameters)
            node.private = event
            rptree.eventDict[node] = event
        end
        node = getDepthFirstNextRight(rptree.treedata, node)
    end
    @printf stdout "\n nbstored = %d \n" length(rptree.eventDict)
end     # end of fillSplittingInfo


"""
This function is used a posteriori to get
    * number of split by diameter and projection
    Presently it computes statistics on leaves diameter
    an return leaves = Array{TreeNode{KeyVector,RPTreeEvent}}() and leafDiameters = Array{Float64,1}()

"""

#    We just scan event dictionary

function analyzeSplittingInfo(rptree::RPTree)
    leafDiameters = Array{Float64,1}()
    leaves = Vector{TreeNode{KeyVector,RPTreeEvent}}()
    nbsplitDiam = 0
    nbsplitProj = 0
    sampleSize = 512
    nbleaf = 0
    #
    if length(rptree.eventDict) == 0
        fillSplittingInfo(rptree)
    end
    #
    res_iter = iterate(rptree.eventDict)
    while res_iter != nothing
        item = res_iter[1]
        #        item[1] is a node, item[2] a RPTreeEvent
        diams  = item[2].diameters
        if item[2].split == splitDiam
            nbsplitDiam = nbsplitDiam+1
        elseif item[2].split == splitProj
            nbsplitProj = nbsplitProj+1
        elseif item[2].split == splitNull
            # we have a leaf
            nbleaf = nbleaf + 1
            leaf = item[1]
            # diameters data have been stored in by fillSplittingInfo
            medDiam = leaf.private.diameters[1]
            #
            push!(leafDiameters, medDiam)
            push!(leaves, leaf)
        end
        res_iter = iterate(rptree.eventDict, res_iter[2])
    end
    meanDiam = mean(leafDiameters)
    @printf stdout "\n meanDiameter at max depth = %f" meanDiam
    @printf stdout "\n nb split diam = %d" nbsplitDiam
    @printf stdout "\n nb split proj = %d\n" nbsplitProj
    # 
    return leaves, leafDiameters
end    # end of analyzeSplittingInfo





"""
    `function fillItemLeafDict(rptree::RPTree)`

    A function to fill a dictionary giving the leaf in which
    is initial item was dipatched
     return a `Dict{Int64,TreeNode}`
    """

function fillItemLeafDict(rptree::RPTree)
    itemLeafDict = Dict{Int64, TreeNode{KeyVector, RPTreeEvent}}()
    # must check that eventDict is filled (or fillSplittingInfo has been called)
    if length(rptree.eventDict) == 0
        fillSplittingInfo(rptree)
    end
    # loop from First left Leaf
    leaf = getFirstLeftLeaf(rptree)
    ileaf = 1
    while leaf != nothing
        # transfer dict from leaf.data to
        idx = collect(keys(leaf.data))
        for i in idx
            itemLeafDict[i] = leaf
        end
        # 
        leaf = getNextLeafRight(leaf)
    end
    return itemLeafDict
end

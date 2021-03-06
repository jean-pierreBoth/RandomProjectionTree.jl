#
#  Simple Tree container with no Key and no Rebalancing
#  To be used in RandomProjectionTree
#
#



using Printf

export
# TreeNode methods
    TreeNode,
    NullNode,
    dumpPos,
    isvoid,
    hasChildren,
    insertChildInNode,
#
#             Tree
#             depth first , left to right iterator
    Tree,
    getFirstLeftLeaf,
    getNextLeafRight,
    getNextLeafLeft,
    getDepthFirstNextRight,
    getDepthFirstNextLeft,
    getLeftSibling,
    getRightSibling

###################################



const debugLevel = 0

"""
# TreeNode{T,U}

FIELDS
------

a Treenode consists in

*  an array of data of generic struct T for which there must exist a Distance.
*  depth : depth in tree , root node is at depth 0.
*  rankInParent : where a node is in parent.children field
*  possibly a parent node
*  an array of children node.
*  some private data of type U for possible use by application using Tree.
    
CONSTRUCTORS
 -----------

1.    `function TreeNode{T,U}(parentArg::TreeNode{T,U}, dataArg::T)`.  
       full constructor with parent node specified and data
    
2.    `function TreeNode{T,U}(dataArg::T, privatedata::U)`. 
       constructor for root node (without parent).

3.    `function TreeNode{T, U}(dataArg::T)`
"""
mutable struct TreeNode{T, U} 
    # fields   
    data::T   # a vector of Data
    depth::Int64
    rankInParent::Int64
    parent::Union{TreeNode{T,U}, Nothing}  # at depth 0 we have a null parent
    # lateral chaining
    children::Array{TreeNode{T,U},1}
    # to store any application dependant data if necessary
    private::Union{U, Nothing}
    #
    l::Threads.ReentrantLock
    # constructor with a parent
    function TreeNode{T,U}(parentArg::TreeNode{T,U}, dataArg::T) where {T,U}
    #    @debug " TreeNode constructor with parent"
        # should check all array have the same dimension
        parent = parentArg
        depth = parent.depth+1
        rankInParent = length(parent.children)+1
        children=Array{TreeNode{T, U},1}()
        nullprivatedata = nothing
        new(dataArg, depth, 0 , parent, children, nullprivatedata, Threads.ReentrantLock() )
        #
        # we could do:
        # n.rankInParent = length(parent.children)+1
        # push!(parent.children, n)
        # n 
        # This last line is necessary so that the last evaluated line is the new object
        # otherwise it is parent.children that is returned !!!
    end
    # constructor for root node of Tree with private data
    function TreeNode{T, U}(dataArg::T, privatedata::U) where{T,U}
    #    @debug " TreeNode constructor root"
        depth = 0
        rankInParent = 0
        parent = nothing
        data=dataArg
        children=Array{TreeNode{T, U},1}()
        new(data, depth,rankInParent, parent, children, privatedata , Threads.ReentrantLock())
    end
    # constructor for root node of Tree without private data
    function TreeNode{T, U}(dataArg::T) where{T,U}
    #    @debug " TreeNode constructor root"
        depth = 0
        rankInParent = 0
        parent = nothing
        data=dataArg
        children=Array{TreeNode{T, U},1}()
        nullprivatedata = nothing
        new(data, depth,rankInParent, parent, children, nullprivatedata, Threads.ReentrantLock())
    end
end


const NullNode = nothing



##############################################################################
#
#  T will most often be an Array of Float64.
#  At least it must be on object for which there is a distance
#
#  All iterations functions take a node as input and return a Nullable
#
############################################################################


"""
 Tree{T,U}  and constructors

a Tree consists in

*  a root node TreeNode{T, U}
*  a lock to serialize operations on Tree if parallel computations are done

"""
mutable struct Tree{T, U}
    #
    root::TreeNode{T, U}
    # to serialize operations on tree
    l::ReentrantLock
    # constructor
    function Tree{T,U}(dataArg::T)  where {T,U}
        root = TreeNode{T,U}(dataArg)
        new(root,ReentrantLock())
    end
end


############################################################################################


"""
dump recursive posiiton

"""
function dumpPos(out::IOStream, n::TreeNode)
    if n.depth > 1
        @printf out "\n depth %d rank in Parent %d" n.depth n.rankInParent
        if n.parent !== nothing
            dumpPos(out, n.parent)
        end
    end
end


function dumpPos(n::TreeNode)
    if n.depth == 0
        @printf stdout "\n\n"
    end
    if n.depth >= 0
        @printf stdout "\n depth %d rank in Parent %d" n.depth n.rankInParent
        if n.parent !== nothing
            dumpPos(n.parent)
        end
    end
end



function isvoid(n::TreeNode)
    if length(n.data) == 0 && length(n.children) == 0
        return true
    else
        return false
    end
end



function hasChildren(n::TreeNode)
    if length(n.hasChildren) == 0
        return false
    else
        return true
    end
end




# This is useful only if we TreeNode constructor do not modify parent
# this is preferred to ensure serialization.
# But that put a lock on all insertion instead of insertion in a given node 

function insertChildInNode(n::TreeNode, child::TreeNode)
    lock(n.l)
#    #
    push!(n.children,child)
    child.parent = n
    child.rankInParent = length(n.children)
    #
    unlock(n.l)
end



function getParent(n::TreeNode)
    return n.parent
end




#
#  find next to the right with first going up until a move right is possible
#  and does the right move.
#  to be abstracted in an iterator Depth first Left To Right
# 

function goToNextRight(t::Tree, n::TreeNode)
    next = nothing
    #
    if debugLevel >= 1
        println("\n\n In goToNextRight :")
        dumpPos(n)        
    end
    #
    node = n
    more = true
    while more 
        if debugLevel >= 2
            dumpPos(node)
        end
        #
        if node.depth != 0 && node.rankInParent < length(node.parent.children)
            more = false   # we can go right
            next = node.parent.children[node.rankInParent+1]
        else
            if node.parent === nothing
                more = false
                next = nothing
            else
                node = node.parent
            end
        end
    end  # end while
    if debugLevel >= 1
        println("\n exiting  goToNextRight :")
        if next !== nothing
            dumpPos(next)
        else
            println("  goToNextRight returns null")
        end
    end
    #
    return next
end


# returns Union{TreeNode{T,U}}, Nothing}
# find next left node by going upward until we can go left.

function goToNextLeft(t::Tree, n::TreeNode)
    next = nothing
    #
    node = n
    more = true
    while more 
        #
        if node.depth != 0 &&  node.rankInParent > 1
            more = false   # we can go left
            next = node.parent.children[node.rankInParent-1]
        else
            if node.parent === nothing
                more = false
                next = nothing
            else
                node = node.parent
            end
        end
    end  # end while
    #
    return next
end


"""
#   function getDepthFirstNextRight(t::Tree, n::TreeNode)

A function to iterate through the tree by depth first left to right traversal

    Returns a node or nothing if end of Tree is reached.
"""
function getDepthFirstNextRight(t::Tree, n::TreeNode)
    next = nothing
    # try to go down, if not possible go right
    if length(n.children) > 0
        next = n.children[1]
    else
        next = goToNextRight(t,n)       # we must go up until some node is not the rightmost in parent
    end
    #
    return next    
end



# iteration depth first right to left

"""
#   function getDepthFirstNextLeft(t::Tree, n::TreeNode)

    A function to iterate through the tree by depth first right to left traversal

    Returns a node or nothing if end of Tree is reached.
"""
function getDepthFirstNextLeft(t::Tree, n::TreeNode)
    next = nothing
    # try to go down, if not possible go right
    if length(n.children) > 0
        next = n.children[length(n.children)]
    else
        next = goToNextLeft(t,n)       # we must go up until some node is not the rightmost in parent
    end
    #
    return next    
end


"""
#   function getLeftSibling(t::Tree, n::TreeNode)
    
    return left sibling (left neighbour at same depth) of a node or nothing
"""
function getLeftSibling(t::Tree, n::TreeNode)
    if n.depth == 0
        return nothing
    end
    leftSibling =  goToNextLeft(t,n)  
    if leftSibling === nothing
        return leftSibling
    end
    #
    # now we must go down to depth of n or until we are at a leaf
    #
    rnode = leftSibling
    if rnode.depth < n.depth
        moreDown = true
        while moreDown
            if length(rnode.children) > 0
                rnode  =  last(rnode.children)
                if rnode.depth == n.depth
                    moreDown = false
                end
            else
                moreDown = false
            end
        end
    end
    #
    if rnode.depth == n.depth
        return rnode
    else
        return nothing
    end
end



"""
#   function getRightSibling(t::Tree, n::TreeNode)

    return right sibling (right neighbour at same depth) of a node or nothing 
"""
function getRightSibling(t::Tree, n::TreeNode)
    if n.depth == 0
        return nothing
    end
    rightSibling = goToNextRight(t,n)  
    if rightSibling === nothing
        return rightSibling
    end
    #
    # now we must go down to depth of n or until we are at a leaf
    #
    rnode = rightSibling
    if rnode.depth < n.depth
        moreDown = true
        while moreDown
            if length(rnode.children) > 0
                rnode = first(rnode.children)
                if rnode.depth == n.depth
                    moreDown = false
                end
            else
                moreDown = false
            end
        end
    end
    #
    if rnode.depth == n.depth
        return rnode
    else
        return nothing
    end
end

       

# to be abstracted in a leaf  left to right iterator

"""
#    function getFirstLeftLeaf(t::Tree)

returns the left most leaf of the tree

to be chained with :
function getNextLeafRight(t::Tree, n::TreeNode)

"""
function getFirstLeftLeaf(t::Tree)
    # iterates depth first until no more children
    node = t.root
    while length(node.children) > 0
        node = node.children[1]
    end
    #
    return node
end

"""
#    function getFirstRightLeaf(t::Tree)
        
returns the rightmost leaf of the tree. 
The call to this functio is to be chained with :  
    function getNextLeafLeft(t::Tree, n::TreeNode)
"""
function getFirstRightLeaf(t::Tree)
    # iterates depth first until no more children
    node = t.root
    while length(node.children) > 0
        node = node.children[length(node.children)]
    end
    #
    return node
end


"""
#  function getNextLeafRight(t::Tree, n::TreeNode)

returns leaf or nothing

"""
function getNextLeafRight(t::Tree, n::TreeNode)
    nextleaf = nothing
    # go up until I can go right
    nextleaf = goToNextRight(t,n)
    # then go down until bottom
    while nextleaf !== nothing && length(nextleaf.children) > 0
        nextleaf = nextleaf.children[1]
    end    
    #
    return nextleaf
end


"""
#  function getNextLeafLeft(t::Tree, n::TreeNode)

returns leaf or nothing
"""
function getNextLeafLeft(t::Tree, n::TreeNode)
    nextleaf = nothing
    # go up until I can go right
    nextleaf = goToNextLeft(t,n)
    # then go down until bottom
    while nextleaf !== nothing && length(nextleaf.children) > 0
        nextleaf = nextleaf.children[length(nextleaf.children)]
    end    
    #
    return nextleaf
end



function deleteSubTree(t::Tree, n::TreeNode)
    parent = n.parent
    nb = length(parent.children)
    if nb > 1
        for i in n.rankInParent+1::nb
            parent.children[i-1] = parent.children[i]
            parent.children[i-1].rankInParent = i-1
        end
        resize!(parent.children, nb-1)
    else
        resize!(parent.children, 0) 
    end
    finalize(n)
end




#!/usr/bin/python -u
import sys, string

###########
## Mo's Revenge programming challenge solution
###########

###########
## Summary:
##   As a code and documentation sample, this file contains more than just the code
##   necessary for solving the problem in an optimal way, using sorted stream vertex and
##   cost data.  The following algorithms and functions are probably the most interesting:
##
##   (1) The algorithm used to determine the next non-overlapping stream in a *sorted*
##       list of streams is O(Mn), where M is the average number of overlapping streams
##       per stream in the list and is independent of n.  Assuming very large n and
##       well-distributed streams, M should be a small constant, reducing O(Mn) to O(n).
##       The stream taken/not taken decisions in (2) below require references to next
##       non-overlapping streams.  Search for function findNextNOStreams().
##
##   (2) The algorithm used to calculate the minimum cost (taken/not taken) decisions
##       through a list of *sorted* stream objects is linear O(n), taking one pass through
##       the list from end to start.  Search for functions calcCostAtEachStream() and
##       calcAccumulatedCostAtStreamIndex().
##
##   (3) The algorithm used to derive the minimum energy is O(1) once the cost minimization
##       decisions at each stream have been computed, and the algorithm used to re-create
##       the optimal path is at most O(n), taking one pass through the stored decisions in
##       tandem with stored references to next non-overlapping streams.  Search for functions
##       deriveOptimalPathFromResults() and deriveOptimalEnergyCostFromResults().
## 
##   As noted in more detail below, while the other costs are approximately O(n),
##   the sorting method chosen to sort the streams in topological order is the built-in
##   Quicksort sort() routine, typically O(n log n).  If the streams are pre-sorted, then
##   it isn't required, potentially saving the O(n log n) step.  It *could* be upgraded
##   to a linear radix sort at the cost of additional complexity.
##       
##   The file contains a preliminary analysis of the problem, some of the thought
##   processes leading to the one-pass solution, lots of code comments, and plenty of
##   checks and bounds for file issues and broken assumptions about streams.  It also
##   includes extra *optional* code for double-checking the one-pass solution with a
##   brute force, graph-based, standard cheapest path algorithm.
###########


# Analyzing the problem:
#
# - It looks like a classical job scheduling problem, except with designated weights on the "jobs" (jet streams)
#   and an added default weight for "no job" (no jet stream).
#
# - The "greedy" algorithm which solves the unweighted scheduling problem by choosing the jet stream which ends
#   first will NOT work for this problem due to the non-uniform stream weights and non-zero default weight.
#
# - A brute force attack of the problem would try to find the shortest path through all possible paths using weighted
#   topological sort (O(n)) (it's a DAG by definition) or Dijkstra's algorithm (basic version is O(n^2)), requiring a
#   very large graph adding edges of default weight between all non-overlapping edges.  Forming a complete graph is at
#   least O(n^2).  Using Dijkstra's algorithm on a complete graph would be very straightforward, but would take
#   *much* longer to evaluate.  At most, it could be used to validate a better, linear algorithm that doesn't
#   require building a complete graph.

# Observations:
#
# - The vertices don't need to be connected for the classical solution in which choosing "jobs" is always
#   desirable.  We can make the edge list match that condition if we discard any jet streams that would cost
#   more for Mo to use than the default.
#
# - We still need to optimize the energy cost over any selection of jet streams while ensuring we choose jet
#   streams as much as possible.
#   - What about searching the permutations of all paths from start to finish?  Unfortunately, determining
#     the best job "x" to place in the optimum job schedule would require recursing through all possibilities
#     at a point in the schedule (at least n calls to start).  No.  It's almost as bad as a complete graph.
#   - Looking at each jet stream, Mo will either use it or not use it.  If he uses it, then we know he cannot
#     use any overlapping jet stream.  So, we could consider throwing out overlapping jet streams when Mo chooses
#     a jet stream and recurse the problem over the remaining jet streams, or throw out the unused jet stream and
#     recurse over the remaining jet streams.  The base case would feed results back through the stack.  We could
#     try both possibilities for each jet stream and choose the best result.
#     Unfortunately, that approach is also bad, since an initial jet stream would have 2*(n-1) recursive calls, with
#     more recursion at each stream, and n decreasing as we progress through the schedule.  Note that the problem
#     decreases in size near the end of the schedule with fewer recursive calls.
#   - The number of possible paths from any stream diminish as the current jet stream (j) moves from 1 to n.
#     Once the streams are sorted topologically, later streams feed back decisions to earlier streams.  We could
#     save time starting at the end point and iterating back to the start point, saving intermediate "relaxation"
#     results on each stream to avoid having to recalculate them.

# Forming the linear algorithm:
#
# - Having decided to "relax" each stream in reverse order and save the results, what do we need to do?
#   - We want to minimize energy cost on the jetstreams, but not using them is worse.  A minimization
#     function must not choose to reject a jetstream because it adds cost, so it must include the default cost.
#     (If we maximize the inverse, 0 is still larger than a negative number, so simply inverting won't work.)
#
#   - At the end of Mo's journey, working backwards, start with accumulated_cost = 0.
#
#   - At each stream we want to minimize the cost of using that stream or not using that stream.  Note that
#     if we use the stream, the next possible stream choice is the next non-overlapping stream.  If we don't
#     use the stream, the next possible stream choice is the next stream in the sorted list.
#
#     Accumulated cost of using the stream (at its start vertex):
#        cost(using stream i)
#      + default_cost(between end of stream i and start of next non-overlapping stream j)
#      + accumulated_cost(at next non-overlapping stream j)
#
#     Accumulated cost of not using the stream (at its start vertex):
#        default_cost(between start of stream i and start of next stream i+1 in list)
#      + accumulated_cost(at next stream i+1 in list) 
#
#     We want to minimize the energy cost, so the accumulated_cost for stream i should be the minimum of
#     the two above functions.  Store the chosen accumulated_cost somewhere, along with the yes/no decision,
#     since we will need to use it (in some way) to construct the optimal path through the streams.
#
#     Since we're working in reverse order, after we complete the base cases, the accumulated_cost results are
#     already calculated and can be retrieved from storage.
#
#     Algorithm complexity to this point:
#       - Sort streams - O(n log n), unless we choose to use radix sort (probably not worth it unless n is *very* large)
#       - Build array of next non-overlapping stream pointers:
#         The complexity is O(Mn), where M is the average number of overlapping streams for each stream.
#         If most to all streams overlap, then worst case is O(n^2), BUT it's very unlikely.
#         Most likely, a stream will overlap with some manageable, random handful of following streams, unrelated to n.
#            When n is large and streams are well-distributed, M will be small in comparison.
#            Typical case for large n will be O(Mn) -> O(n).
#       - Compare taking a stream with not taking a stream and store results - O(n).
#     Worst shown, assuming large n and well-distributed streams:  O(n log n) for the initial sort, unless we
#     use radix sort (probably not worth the extra complexity for O(pk + pn), as O(n log n) optimized and pre-compiled
#     Quicksort is already built-in to Python).

# Deriving total energy cost:
#
# - At the end of the algorithm, the total energy cost at the first taken stream's start vertex should be available.
#   If the start vertex is not zero, then the final cost should include the default cost from point zero to the
#   start of the first taken stream.

# Deriving the optimal path:
#
# - We cannot just look at the yes/no results for each stream to derive the optimal path.  In most cases, the local
#   decision to take or not to take the stream will be -yes-.
#   - Looking at how more complicated graph algorithms work, one constructs the desired path by following optimal
#     parent pointers through the graph.  In this case, we don't build a connected graph, but we know we usually want
#     to take the next available stream.  If we take a stream, that's the next non-overlapping stream; otherwise, it's
#     the next stream in the ordered list.
#   - The first stream with the final accumulated cost (excluding any preceding default cost) should indicate if
#     it was chosen as taken or not taken.  If taken, look at the decision for the next non-overlapping stream;
#     if not taken, look at the decision for the next stream in the ordered list.  O(n) complexity.
#   - If the first stream is taken in both sample files, additional confidence that the first stream will have
#     valid taken/not taken results can be supplied by adding a very short and minimally cost-effective stream
#     beginning just before and barely overlapping the first stream.

# Additional OPTIONAL validation:
#
# - If we build a complete, explicit graph of all possible paths, we can run Dijkstra's well known shortest-path
#   algorithm on the graph.  We can traverse the parent pointers to re-create the path and look at the accumulated
#   energy cost at the final vertex.  Building the entire graph is straightforward but O(n^2), so it's very slow
#   with large n.  Dijkstra's basic algorithm is also O(n^2); it could be improved using Fibonacci heaps to
#   O(m + n log n), but the added complexity would also need to be validated.  The same applies to any simplifications
#   of the algorithm due to assumptions (e.g., DAG).
# - The optional validation can be waived given enough confidence in the "one-pass" solution.


#####################
#### CODE BEGINS ####
#####################

# Global options to set through command line arguments.  They can be hardwired here, but it's not necessary.
# Do NOT choose to validate with Dijkstra unless you want to wait a very long time for large n.
do_debug = False                       # print extra debugging info
do_val_with_Dijkstra = False           # also build explicit DAG - O(n^2)!!! - and validate results with Dijkstra


# MoStream class
#   It holds and accesses the start, end, cost, accumulated cost, next non-overlapping stream pointer,
#   and yes/no decisions for a stream.
#   When instantiating the class, provide the start, end, and cost for the stream.
#   End is assumed to be >= start, and cost is assumed to be >= 0.
#
#   Usage:  myNewStream = MoStream(start_vertex, end_vertex, total_stream_cost)
class MoStream:
    # global attributes for the class
    start = 0            # start mile marker
    end   = 0            # end mile marker
    cost  = 0            # cost over entire stream
    accumCost = 0        # accumulated cost to be calculated as we progress
    nextNOStream = None  # reference to next non-overlapping stream (as an index into an array of streams)
    useStream = False    # flag to indicate resolution of min cost function

    # Initialize a new MoStream object, given start vertex, end vertex, and stream cost values
    def __init__(self, inputStart, inputEnd, inputCost):
        self.start = inputStart    # set start mile marker
        self.end = inputEnd        # set end mile marker
        self.cost = inputCost      # set cost over entire stream
        if(do_debug):
            print 'DEBUG: Stream initialized with start %i, end %i, cost %i' % (self.start, self.end, self.cost)

    # Boolean - Check if another stream overlaps with this one
    # Usage:  streamsOverlap = stream.doesStreamOverlap(other_stream)
    #   Is the other stream's start between this stream's start and end, or
    #   is the other stream's end between this stream's start and end?
    #   - We don't care if the vertices are equal at start or end; we assume Mo can transition cleanly.
    def doesStreamOverlap(self, otherStream):
        if(((otherStream.start > self.start) and (otherStream.start < self.end)) or 
           ((otherStream.end > self.start) and (otherStream.end < self.end))):
            return True
        return False

### END OF MoStream CLASS ###


################
##### EXTRA CLASSES FOR VALIDATION USE ONLY - MoGraphEdge, MoGraph  
################

### For OPTIONAL validation ONLY with Dijkstra ###
### This is not used with the simplified "one-pass" algorithm.
# MoGraphEdge class
#   It holds and accesses the start, end, cost, next edge pointer, and default (is/is not a stream)
#   values for an edge in the connected DAG, used for validation.
#   End is assumed to be >= start, and cost is assumed to be >= 0.
#
#  Usage:  myNewEdge = MoGraphEdge(start_vertex, end_vertex, stream_cost, is_default_path)
#  For streams, set is_default_path = 0.  For edges connecting streams, set is_default_path = 1.
class MoGraphEdge:
    start = 0            # start mile marker
    end   = 0            # end mile marker
    cost  = 0            # cost over entire edge
    nextEdge = None      # pointer to next MoGraphEdge in a linked list
    default = True       # is this edge connecting streams (default) or a stream itself (not default)

    #  For streams, set isDefault = 0.  For edges connecting streams, set isDefault = 1.
    #  Cost should be the total cost for the stream (for default, calculate it over the entire distance).
    def __init__(self, inputStart, inputEnd, inputCost, isDefault):
        self.start = inputStart
        self.end = inputEnd
        self.cost = inputCost
        self.default = isDefault
#        if(do_debug):
#            print 'DEBUG: Edge initialized with start %d, end %d, cost %d, default %i' % (self.start, self.end, self.cost, self.default)

### END OF MoGraphEdge CLASS ###


### For OPTIONAL validation ONLY with Dijkstra ###
### This is not used with the simplified "one-pass" algorithm.
# MoGraph class
#   It holds and accesses the graph's edges, the outbound connection count for each vertex, a total
#   count of vertices and edges, an optional "directed" flag, the default cost per unit distance
#   for edges connecting streams, and a lookup dictionary for start vertices.
#
#  Usage:  myNewGraph = MoGraph(default_non_stream_cost)
class MoGraph:
    graphEdges = []          # Array of linked lists of edges connected to each vertex; indexed by vertex
    nodeDegrees = []         # Array of counts of outgoing edges for each vertex; indexed by vertex
    numVertices = 0          # Total number of vertices in the graph
    numEdges = 0             # Total number of edges in the graph (yes, it will be O(n^2))
    directed = True          # It's a DAG by definition
    defaultCost = 0          # Holder variable for default unit travel cost
    nodeLookupDict = {}      # Hash for mapping start vertex values to array indices used by graphEdges and nodeDegrees

    # Initialize new MoGraph object with given default cost per unit travel
    def __init__(self, inputDefaultCost):
        self.defaultCost = inputDefaultCost

    # Number value - Calculate the total default cost over a non-stream travel segment,
    #                given the starting and ending vertices of the segment.
    # Usage:  local_cost = graph.getDefaultCostForDistance(start, end)
    def getDefaultCostForDistance(self, start, end):
        totalCost = (end - start) * self.defaultCost    # cost = (distance) * (unit_cost)
        if(totalCost < 0):                              # If end < start, we'll need to invert cost
            return -(totalCost)
        return totalCost
    
    # Void - Add an edge to the graph, given start and end vertices, the cost for the edge,
    #        and a boolean indicator of stream/default.
    # NOTE:  When isDefault=True, the cost will be calculated by the getDefaultCostForDistance
    #        function and the cost argument will be ignored.
    def addEdge(self, start, end, cost, isDefault):
        if(self.nodeLookupDict.has_key(start)):          # If we've already seen the start vertex,
            nodeIndex = self.nodeLookupDict[start]       # look up its vertex index mapping
        else:
            nodeIndex = self.numVertices                 # It's new; the new last index will be the index count
            self.nodeLookupDict[start] = nodeIndex       # Add the mapping
            self.nodeDegrees.append(0)                   # Add the node degrees value (starts at 0 outgoing connections)
            self.graphEdges.append(None)                 # Add an empty linked list of connected edges
            self.numVertices = self.numVertices + 1      # Increment the number of total vertices
        if(isDefault):                                              # If an edge connecting streams,
            edgeCost = self.getDefaultCostForDistance(start, end)   # calculate the default cost over the length
        else:
            edgeCost = cost                                         # Otherwise, use the given cost for the stream
        newEdge = MoGraphEdge(start, end, edgeCost, isDefault)      # Create the new MoGraphEdge object
        headEdge = self.graphEdges[nodeIndex]                       # Get the head of the linked list of edges for this vertex
        if(headEdge != None):                                       # If one exists, link it to the new edge
            newEdge.nextEdge = headEdge                             
        self.graphEdges[nodeIndex] = newEdge                        # Update the head of the linked list to the new edge
        self.numEdges = self.numEdges + 1                           # Increment the number of total edges
        self.nodeDegrees[nodeIndex] = self.nodeDegrees[nodeIndex] + 1  # Increment the number of outgoing edges for this vertex

    # Void - Add a final vertex to the graph with Mo's end location if it doesn't already exist.
    #        As defined by the problem, the end vertex should have no outgoing edges, and it will
    #        not likely be added by the iterative graph building algorithm.
    def addFinalVertex(self, end):
        if(not self.nodeLookupDict.has_key(end)):
            nodeIndex = self.numVertices
            self.nodeLookupDict[end] = nodeIndex
            self.nodeDegrees.append(0)
            self.graphEdges.append(None)
            self.numVertices = self.numVertices + 1

### END OF MoGraph CLASS ###


################
##### MAIN CLASS FOR PROCESSING STREAMS AND PRINTING RESULTS - MoRev Class
################

# MoRev class
#   This class parses the input file, creates the streams, runs the algo, and prints the results as directed
#   by the problem statement.  If the optional validation flag has been set, it will also run the very long
#   (but straightforward) shortest path algorithm and compare the results to the one-pass algorithm.
#
#   When instantiating the class, provide the input file from the command line.
#   Usage:  moRevInstance = MoRev(input_file_pathname)
class MoRev:
    streamFile = ''                  # store the input file string
    defaultCostPerMile = 0           # store the default cost per mile, held in the input file
    streamDict = {}                  # this dictionary will hold all the valid streams
    streamList = []                  # this list contains sorted streams in ascending mile marker order (sorted eventually)
    listRequiresSort = True          # indicate whether the stream list is pre-sorted or not during parsing
    lastMileMarker = 0               # this holds the last mile marker (greatest end value for any stream)
    linearTuples = []                # path tuple list for one-pass algorithm
    linearEnergy = 0                 # total energy cost for one-pass algorithm
    # Instance variables used only for validation
    graph = None                     # MoGraph object - USED ONLY WHEN VALIDATION IS ENABLED
    dijkstraTuples = []              # path tuple list for validation algorithm - USED ONLY WHEN VALIDATION IS ENABLED
    dijkstraEnergy = 0               # total energy cost for validation algorithm - USED ONLY WHEN VALIDATION IS ENABLED


    # Returns new MoRev instance, parses the input file, runs the one-pass algorithm and the optional
    # validation algorithm with checks, and prints the results.
    def __init__(self, inputFile):
        self.streamFile = inputFile           # store the input file string
        if(not (self.streamFile == '')):      # parse file if a name has been supplied, or exit with error
            success = self.parseStreamFile()
        else:
            print 'You must provide an input file to this program.  Usage:  python mo_rev.py <input_file_name>'
            exit(1)
        if(not success):                      # exit with error if any problem opening or using the file
            print 'The parser failed to process the stream input file.  Exiting...'
            exit(1)

        if(self.listRequiresSort):    # if the streams in streamList have not been pre-sorted,
            self.sortStreamList()     #   sort the list of keys into the dictionary by start points
                                      #   (and ends, if starts are equal) - parser checks order; may skip O(n log n)
        self.findNextNOStreams()      # fill in the nextNOStream element of each stream,
                                      #   indicating the next following non-overlapping stream's index into the stream list
        self.calcCostAtEachStream()   # calculate the cost at each stream
        self.printChosenStreamsAndCost()  # print the total cost and stream tuples for the optimal path

        # If the optional validation flag has been set, run the VERY long O(n^2) graph build process and
        # run Dijkstra's algorithm on it.  Compare its results to the one-pass algorithm's results above.
        if(do_val_with_Dijkstra):
            self.buildGraph()
            self.runDijkstraOnGraph(self.graph)
            self.compareLinearToDijkstra()


    # Boolean - Parse the supplied stream file, setting default cost and building an initial, un-sorted stream list.
    #           Returns True if the function succeeded; returns false for file access errors or badly formatted data.
    #           It expects self.streamList to be an empty list.  Due to the number of checks on the parsed data,
    #           parsing has been split into multiple functions in a couple of layers under parseAndAddStreamsFromLines.
    def parseStreamFile(self):
        if(not self.streamFile):       # can't open a non-existent file
            return False
        if(len(self.streamList) > 0):  # stream list should be empty on object initialization
            print 'Parser error:  The parser should start with an empty stream list.  Has it been called multiple times?'
            return False
        try:
            fileInput = open(self.streamFile, 'r')   # open the stream file for reading
        except:                                      # if we didn't open the file, then exit and signal error
            print 'Parser error:  Could not open the supplied file %s.' % (self.streamFile)
            return False

        allLines = fileInput.readlines()    # read all the lines into separate strings; we're assuming the file isn't enormous
        fileInput.close()                   # close the file
        parseStatus = self.parseDefaultCostFromLine(allLines[0])   # default cost should be in the first line (0)
        if(not parseStatus):
            return False
        return self.parseAndAddStreamsFromLines(allLines, 1)       # streams should be in remaining lines (from 1)


    # Boolean - This function takes a string expected to contain a single numerical value,
    #           the default cost per unit distance when not using a stream.  It returns
    #           False if the string or the default cost does not match expectations.
    #           Otherwise, it sets the defaultCostPerMile variable and returns True.
    def parseDefaultCostFromLine(self, lineString):
        columns = string.split(lineString)  # get the separated values from the first line
        numColumns = len(columns)           # count the columns

        # The first line should have only one number in it
        if(not(len(columns) == 1)):
            print 'Parser error:  The first line does not contain a single value for the default cost per mile.'
            return False  
        costPerMile = eval(columns[0])                 # default cost per mile is the only value in the line
        if((not costPerMile) or (costPerMile == 0)):   # if it doesn't exist or make sense, exit with an error
            print 'Parser error:  The default cost per mile is undefined or zero.'
            return False
        self.defaultCostPerMile = costPerMile          # save value in instance variable if valid and return true
        return True


    # Boolean - This function takes the array of lines from the stream input file and
    #           an index pointing to the first line of the array with stream inputs, and
    #           it returns false if it encounters any errors in the data while parsing.
    #           It parses the start, end, and cost values from each line, then calls other
    #           functions to check the values and to add them to the instance's streamList variable.
    def parseAndAddStreamsFromLines(self, allLines, startLine):
        # Determine if the stream list needs to be sorted (save a O(n log n) sort)
        localRequiresSort = False           # start with assumption that the list is pre-sorted
        (lastStart, lastEnd) = (-1, -1)     # keep track of previous values to determine if sorted

        # Parse the start, end, and cost values from all lines from startLine
        for line in allLines[startLine:]:   # cycle over the rest of the lines with stream definitions
            columns = string.split(line)    # split the line into separate numbers
            numColumns = len(columns)       # check column count
            if(numColumns > 0):             # we'll skip a line if it's empty, assuming that might happen
                if(numColumns < 3):         # if fewer than 3 numbers, the line is invalid, so exit with an error
                    print 'Parser error:  Each line of the input file must contain 3 values for the stream start, finish, and cost.'
                    return False
                streamStart = eval(columns[0])          # start is first number
                streamEnd = eval(columns[1])            # end is second number
                streamCost = eval(columns[2])           # cost is third number
                if(self.streamValuesHaveErrors(streamStart, streamEnd, streamCost)):    # check for start > end, negative cost
                    return False
                if(self.didAddStream(streamStart, streamEnd, streamCost)):
                    if(streamEnd > self.lastMileMarker):   # update lastMileMarker if end point for this stream is later
                        self.lastMileMarker = streamEnd
                    if(not localRequiresSort):         # only continue to check order if list still doesn't need to be sorted
                        if((streamStart < lastStart) or ((streamStart == lastStart) and (streamEnd < lastEnd))):
                            # If this stream's start is before the last stream's start, or they're equal and
                            # the stream's end is before than the last stream's end, then the list isn't sorted.
                            localRequiresSort = True                              
                            if(do_debug):
                                print 'DEBUG:  Parser determined streams needed to be sorted with sequential streams %d %d and %d %d.' % (lastStart, lastEnd, streamStart, streamEnd)
                        else:                                                     # otherwise, still going to check, so save values
                            lastStart = streamStart
                            lastEnd = streamEnd                    
        # don't update pessimistic instance variable unless we successfully parse all the streams
        self.listRequiresSort = localRequiresSort
        return True


    # Boolean - This function takes a stream start, end, and cost and returns True if the stream
    #           was added to the instance's streamList variable.
    #           The function will ignore streams with the following conditions:
    #           (a) If the stream start and end vertices are identical, then the stream
    #               is zero-distance and pointless.
    #           (b) If the average cost per unit distance of a stream is greater than
    #               the default, then Mo would never choose it.
    #           (c) If a pair of start and end vertices match an already saved pair,
    #               the function will ignore the stream if the new cost is higher.
    #               ** Otherwise, it will overwrite the saved cost with the new, lower cost.
    def didAddStream(self, streamStart, streamEnd, streamCost):
        # Discard any pointless streams of zero length, especially before dividing by length below...
        if(streamStart == streamEnd):
            print 'Parser warning:  Discarding zero-length stream:  start %d, end %d, cost %d...' % (streamStart,
                                                                                                     streamEnd,
                                                                                                     streamCost)
            return False

        # Check cost per mile; discard if above default cost
        costPerMile = streamCost / (streamEnd - streamStart)  # check average cost per mile (integer okay)
        if(do_debug):
            print 'DEBUG: cost per mile %d for stream %d %d %d' % (costPerMile, streamStart, streamEnd, streamCost)

        # If average cost per mile is greater than or equal to the default, throw out the edge (integer okay).
        # We definitely won't use it if it costs more energy to use than the default.
        if(costPerMile >= self.defaultCostPerMile):         
            print 'Parser warning:  Discarding stream with start %d, end %d, cost %d, due to higher cost per mile %d than the default %d...' % (streamStart, streamEnd, streamCost,
                                                                                                                                costPerMile, self.defaultCostPerMile)
            return False
            
        # Otherwise, make a unique key with start and end values appended together to check for duplicates
        streamKey = '%dto%d' % (streamStart, streamEnd)
        stream = MoStream(streamStart, streamEnd, streamCost)
        if(self.streamDict.has_key(streamKey)):
            previousCost = self.streamDict[streamKey].cost
            # If we have a duplicate and its cost is lower, use the lower cost.
            if(streamCost < previousCost):
                print 'Parser warning:  Duplicate start and end vertices %d and %d.  Replacing cost with lower cost %d...' % (streamStart,
                                                                                                                           streamEnd,
                                                                                                                           streamCost)
                self.streamDict[streamKey].cost = streamCost
            else:  # Otherwise, throw out the duplicate with a larger cost (it wouldn't be used).
                print 'Parser warning:  Duplicate start and end vertices %d and %d.  Cost is higher or same.  Ignoring...' % (streamStart,
                                                                                                                           streamEnd)
        else:
            self.streamDict[streamKey] = stream   # add new stream object to dictionary under key
            self.streamList.append(stream)        # add new stream object to end of stream list
            if(do_debug):
                print 'DEBUG: added stream to dictionary and list with key ' + streamKey
            return True
        return False


    # Boolean - This function takes a stream start, end, and cost and returns True if any of the
    #           values violate assumptions about the data.  Currently, the end vertex must be
    #           greater than or equal to the start vertex, and the cost of the stream must be
    #           greater than or equal to zero.
    def streamValuesHaveErrors(self, streamStart, streamEnd, streamCost):
        if(streamEnd < streamStart):
            print 'Parser error:  Assumption that start < end has been violated for stream %d %d cost %d' % (streamStart,
                                                                                                             streamEnd,
                                                                                                             streamCost)                  
            return True
        if(streamCost < 0):
            print 'Parser error:  Assumption that all costs are positive has been violated for stream %d %d cost %d' % (streamStart,
                                                                                                                        streamEnd,
                                                                                                                        streamCost)
            return True
        return False
        
    
    # Void - This sorts the instance's list of stream objects by ascending start point (and end
    #        points, if starts are equal).
    #        It currently uses the built-in Quicksort sort() routine - typically O(n log n).
    #        If we typically encountered excessively large number of streams, we could upgrade
    #        it to linear O(pk + pn) with radix sort, but the proven built-in O(n log n) sort
    #        should be suitable.
    def sortStreamList(self):
        def streamCmpFunction(a, b):
            startCompare = cmp(a.start, b.start)  # Compare start points
            if(startCompare == 0):                # If equal, return comparison of end points
                return cmp(a.end, b.end)
            else:
                return startCompare               # Return start comparison if inequal
        # Quicksort streamList in place using the above function
        self.streamList.sort(streamCmpFunction)
        

    # Void - This function sets the nextNOStream reference in each stream to the streamList index for the
    #        next following stream which does not overlap for each stream in the instance's streamList.
    # NOTE - The streamList[] array must be sorted before calling this function!
    #        The simplified overlap check looks only at the this.end and next.start since
    #        next.start should be >= this.start, and next.end should be >= next.start.
    #        (We could use the more general overlap function inside the MoStream object...)
    def findNextNOStreams(self):
        for streamIndex in range(0, len(self.streamList)):              # over all items in streamList do...
            currentStream = self.streamList[streamIndex]                # get stream at index
            foundNOStream = False                                       # haven't yet found a non-overlapping stream
            nextStreamIndex = streamIndex + 1                           # next stream in series is current + 1
            # While we haven't found a stream and the next stream is still within the bounds of the list,
            # check if the next stream's start is equal to or greater than the end of the current one (non-overlapping).
            # If not, increment next stream and keep checking.
            while ((not foundNOStream) and (nextStreamIndex < len(self.streamList))):
                nextStream = self.streamList[nextStreamIndex]
                if(nextStream.start >= currentStream.end):
                    foundNOStream = True
                else:
                    nextStreamIndex = nextStreamIndex + 1
            # If we found a subsequent non-overlapping stream, set the internal reference to it, else set it to None.
            # We can check for None when calculating the cost.
            if(foundNOStream):
                currentStream.nextNOStream = nextStreamIndex  # save the index to the next non-overlapping stream in this stream
            else:
                currentStream.nextNOStream = None             # there is no next non-overlapping stream

            if(do_debug):
                if(foundNOStream):
                    print 'DEBUG:  Next NO stream for stream index %i is stream index %i' % (streamIndex, nextStreamIndex)
                else:
                    print 'DEBUG:  No next NO stream was found for stream index %i' % (streamIndex)


    # Void - This function iteratively calculates the accumulated cost over the streams in the
    #        instance's streamList, starting from the end.  It assumes the stream list is sorted and
    #        all streams have a reference to the index for the next non-overlapping stream in the list.
    def calcCostAtEachStream(self):
        # end of last stream is the defined end of the journey
        numStreams = len(self.streamList)           
        streamIndex = len(self.streamList) - 1       # start with last stream
        while (streamIndex >= 0):                    # over streams n-1 to 0:
            # get taken/not taken result and accumulated cost for this stream
            (streamChosen, accumulatedCost) = self.calcAccumulatedCostAtStreamIndex(streamIndex)
            self.streamList[streamIndex].useStream = streamChosen        # store taken/not taken
            self.streamList[streamIndex].accumCost = accumulatedCost     # store accumulated cost
            streamIndex = streamIndex - 1                                # decrement stream counter


    # Number value - Returns total cost over an interval between two streams, given start and end vertices
    def defaultCostForInterval(self, start, end):
        cost = (self.defaultCostPerMile * (end - start))
        if(cost < 0):       # really, it should be >= 0, but let's be careful anyway
            cost = -(cost)
        return cost


    # Tuple - (stream_taken, accumulated_cost)
    #         This function assumes the stream list is sorted and all streams have a reference to
    #         the index for the next non-overlapping stream in the list.  It also assumes that
    #         results for streams with higher indices have already been calculated.
    #         This function calculates the minimum accumulated cost for taking or not taking the
    #         stream at streamList[streamIndex].  It returns the boolean decision and the new
    #         accumulated cost at the stream's start vertex.
    def calcAccumulatedCostAtStreamIndex(self, streamIndex):
        # minimize cost of using a stream or not using a stream
        # Accumulated cost of using a stream:
        #   cost(stream n) + defaultCost(between n.end and nextNOStream.start) + accumCost(nextNOStream)
        # Accumulated cost of not using a stream:
        #   defaultCost(between n.start and (n+1).start) + accumCost(n+1))
        currentStream = self.streamList[streamIndex]       # fetch stream at streamIndex
        if(not (currentStream.nextNOStream == None)):      # fetch next non-overlapping stream, if it exists
            nextNOStream = self.streamList[currentStream.nextNOStream]
        else:
            nextNOStream = None
        if(not (nextNOStream == None)):                    # if next NO stream exists, use it in the cost calculation
            costUsingStream = currentStream.cost + self.defaultCostForInterval(currentStream.end, nextNOStream.start) + nextNOStream.accumCost
        else:                                              # otherwise, calculate default cost to last mile marker (the end)
            costUsingStream = currentStream.cost + self.defaultCostForInterval(currentStream.end, self.lastMileMarker)

        if(streamIndex < (len(self.streamList) - 1)):      # if this isn't the last stream, calc not taken cost to next stream
            nextStreamInList = self.streamList[streamIndex + 1]
            costNotUsingStream = self.defaultCostForInterval(currentStream.start, nextStreamInList.start) + nextStreamInList.accumCost
        else:                                              # otherwise, calc not taken cost to the end (should be 0 by definition)
            costNotUsingStream = self.defaultCostForInterval(currentStream.start, self.lastMileMarker)

        if(do_debug):
            print 'DEBUG:  At stream %i, costUsingStream is %d, costNotUsingStream is %d' % (streamIndex, costUsingStream, costNotUsingStream)
        # return result tuple (stream_taken, accumulated_cost)
        if(costUsingStream <= costNotUsingStream):
            return (True, costUsingStream)
        else:
            return (False, costNotUsingStream)

    
    # Tuple - (list_of_chosen_stream_tuples, list_of_indices_to_chosen_streams)
    #         This function assumes the results of the algorithm have been calculated and the individual min cost
    #         taken/not taken decisions reside in each stream's useStream variable.
    #         This function returns a tuple with two list objects.  The first is a list of tuples of start and end
    #         vertices for taken streams (for printing the desired answer).  The second is a list of indices of
    #         taken streams in the optimal path (for ease of verification).
    def deriveOptimalPathFromResults(self):
        chosenTupleList = []              # List of (start, end) stream tuples for optimal path
        chosenIndexList = []              # List of chosen stream indices
        numStreams = len(self.streamList)

        nextNOStreamIndex = 0             # Start with final results stored in the first stream in the list
        while((nextNOStreamIndex != None) and (nextNOStreamIndex < numStreams)):
            nextNOStream = self.streamList[nextNOStreamIndex]
            if(nextNOStream.useStream):              # Add stream info to lists if taken, get next NO stream
                chosenTupleList.append((nextNOStream.start, nextNOStream.end))
                chosenIndexList.append(nextNOStreamIndex)
                nextNOStreamIndex = nextNOStream.nextNOStream
            else:                                    # Otherwise, check the following stream in the list (it should also be NO)
                nextNOStreamIndex = nextNOStreamIndex + 1
                if(do_debug):
                    print 'DEBUG:  Path info: Min cost decision is NOT TAKEN for next non-overlapping stream %d %d %d, cost %d' % (nextNOStream.start, 
                                                                                                                                   nextNOStream.end, 
                                                                                                                                   nextNOStream.cost,
                                                                                                                                   nextNOStream.accumCost)
                    print 'DEBUG:  Moving on to the next stream in the list.'
        return (chosenTupleList, chosenIndexList)


    # Number value - Returns total energy cost as calculated by the algorithm.
    #                This function assumes the results of the algorithm have been calculated and the accumulated
    #                energy cost from the start of the first taken stream to the last mile marker resides in the
    #                first chosen stream's accumCost variable.
    def deriveOptimalEnergyCostFromResults(self, firstUsedStreamIndex):
        firstUsedStream = self.streamList[firstUsedStreamIndex]   # get first taken stream
        energyCost = firstUsedStream.accumCost                    # get preliminary energy cost
        startOfFirstUsedStream = firstUsedStream.start
        if(startOfFirstUsedStream > 0):   # if start vertex is not 0, add the default cost from 0 to the start vertex
            energyCost = energyCost + (self.defaultCostForInterval(0, startOfFirstUsedStream))
        return energyCost


    # Void - This function assumes the results of the algorithm have all been calculated.  It calls
    #        the functions to derive the optimal path and the optimal energy cost through the path.
    #        It checks for errors and then prints the optimal energy cost and tuples with the start
    #        and end vertices for all streams in the optimal path.
    #        NOTE:  It also saves the "linear" one-pass algorithm's results in instance variables
    #               in case full verification with an explicit graph has been requested.
    def printChosenStreamsAndCost(self):
        (chosenTupleList, chosenIndexList) = self.deriveOptimalPathFromResults()   # get optimal path from results
        if((len(chosenTupleList) == 0) or (len(chosenIndexList) == 0)):            # exit if results make no sense
            print 'ERROR:  List of chosen streams for the path is empty!  Exiting...'
            exit(1)

        energyCost = self.deriveOptimalEnergyCostFromResults(chosenIndexList[0])   # get and print optimal energy cost from results
        print 'Minimum total energy required to fly the complete path is %i' % (energyCost) 
        validatedEnergyCost = self.checkCostThroughChosenList(chosenIndexList)     # check energy cost by summing through path
        if(do_debug):
            print 'DEBUG:  The double-checked cost through the path is %i' % (validatedEnergyCost)
        if(energyCost != validatedEnergyCost):                                     # if costs are not equal, signal error
            print 'ERROR:  The accumulated energy cost %d does not match the energy cost %d summed through the chosen path!  Exiting...' % (energyCost, validatedEnergyCost)
            
        print 'The optimal sequence of jet streams is the following:'              # print the optimal path stream vertices
        print chosenTupleList
        print
        # Save references to the one-pass algorithm results for use in OPTIONAL validation
        self.linearTuples = chosenTupleList
        self.linearEnergy = energyCost


    # Number value - The function returns the accumulated cost calculated by summing the costs along the optimal path.
    #                It assumes the results of the algorithm have all been calculated and receives the list of
    #                chosen stream indices.  It steps through the streams and sums the cost of the streams and the
    #                default costs between streams.
    def checkCostThroughChosenList(self, chosenList):
        lastEnd = 0         # start at location 0
        accumCost = 0       # total cost starts at 0
        for streamIndex in chosenList:     # over all chosen streams, add stream cost and cost from lastEnd to current stream's start
            currentStream = self.streamList[streamIndex]
            accumCost = accumCost + self.defaultCostForInterval(lastEnd, currentStream.start) + currentStream.cost 
            if(do_debug):
                print 'DEBUG:  Cost check: index %i start %d last_end %d stream_cost %d, accum_cost %d' % (streamIndex, currentStream.start,
                                                                                                           lastEnd, currentStream.cost, accumCost)
            lastEnd = currentStream.end
        if(do_debug):
            print 'DEBUG:  Cost check: last mile marker is %d' % (self.lastMileMarker)
        accumCost = accumCost + self.defaultCostForInterval(lastEnd, self.lastMileMarker)
        return accumCost


    # Void - FOR VALIDATION PURPOSES ONLY
    #        This function builds a fully connected, explicit, directed graph from the stream list.
    #        The number of edges will be O(n^2), where n is the number of streams.  Use the graph
    #        only for validating results obtained by the one-pass algorithm.
    #        The list of streams should already be pruned for any streams with duplicate vertices
    #        (use one with min cost) and for any streams with costs greater than default over their
    #        lengths (discard as wouldn't be used).
    def buildGraph(self):
        # self.lastMileMarker holds max end point (defined end of graph)
        self.graph = MoGraph(self.defaultCostPerMile)
        # Build edges from 0 to start of every stream (they're all possible, even if most are very unlikely)
        # Also build edges for each stream
        # Also build edges from end of each stream to end of graph (they're all possible, even if most are very unlikely)
        self.graph.addEdge(0, self.lastMileMarker, 0, True)                   # add start to end
        for streamIndex in range(0, len(self.streamList)):                    # over all streams:
            stream = self.streamList[streamIndex]
            # don't want edge from 0 to 0
            if(stream.start > 0):
                self.graph.addEdge(0, stream.start, 0, True)                  # add start to stream.start
            self.graph.addEdge(stream.start, stream.end, stream.cost, False)  # add stream
            # don't want edge from end to end
            if(stream.end < self.lastMileMarker):
                self.graph.addEdge(stream.end, self.lastMileMarker, 0, True)  # add stream.end to end
        # Build edges from end of every stream to start of all following non-overlapping streams
        for streamIndex in range(0, len(self.streamList) - 1):                      # over all (i) but last stream:
            firstStream = self.streamList[streamIndex]
            for secondStreamIndex in range(streamIndex + 1, len(self.streamList)):  # and over following stream (j) to last stream:
                secondStream = self.streamList[secondStreamIndex]
                if(not firstStream.doesStreamOverlap(secondStream)):                # if streams don't overlap, add edge between them
                    # don't want to create edge when nodes match
                    if(secondStream.start > firstStream.end):
                        self.graph.addEdge(firstStream.end, secondStream.start, 0, True)
        # Add last vertex with no directed edges
        self.graph.addFinalVertex(self.lastMileMarker)


    # Void - FOR VALIDATION PURPOSES ONLY
    #        This function takes a fully connected, explicit, directed MoGraph object and runs Dijkstra's
    #        algorithm on it to find the optimal path.  It relaxes the edges by choosing edges with minimum
    #        energy cost.  Distances are included in the costs (weights) and are not separated from them in
    #        the cost function.
    #        NOTE:  The algorithm is the basic, published algorithm without any optimizations in order to
    #               increase validation confidence.
    def runDijkstraOnGraph(self, graph):
        start = 0                  # start with vertex 0
        costs = []                 # array of costs at each vertex
        inTree = []                # array tracking vertex visitation
        parents = []               # array with arrays of parent vertex, edge start vertex, and edge object pointer for each vertex
                                   # - usually parent vertex should be sufficient, but the extra info is used
                                   #   for optionally printing debug information
        largeValue = (1 << 30) - 1 # Dijkstra needs to start search with a very large number; assumes the stream file values are smaller

        # Initialize array tracking list, cost list (to large value), and parents lists
        for currentVertex in range(0, graph.numVertices):
            inTree.append(False)
            costs.append(largeValue)
            parents.append([-1, -1, None])

        costs[start] = 0
        currentVertex = start
        while(inTree[currentVertex] == False):   # if we haven't fully visited this vertex, work on its edges
            inTree[currentVertex] = True
            currentEdge = graph.graphEdges[currentVertex]
            while(currentEdge != None):          # step through all edges leaving this vertex
                candidateNextVertex = graph.nodeLookupDict[currentEdge.end]
                currentCost = currentEdge.cost
                # if cost at next candidate is less than accumulated cost + cost of this edge, use the lower cost
                # and save (or update) the parent information
                if(costs[candidateNextVertex] > costs[currentVertex] + currentCost):
                    costs[candidateNextVertex] = costs[currentVertex] + currentCost
                    parents[candidateNextVertex] = [currentVertex, currentEdge.start, currentEdge]
                currentEdge = currentEdge.nextEdge
            currentVertex = 0                    # start again at vertex 0
            currentCost = largeValue             # start with large base cost
            # over all vertices, for those that haven't been thoroughly visited and have a cost
            # less than the current minimum cost, update the minimum cost and an index to the
            # vertex (with the new minimum) to use on the next loop.
            for vertexIndex in range(0, graph.numVertices):
                if((inTree[vertexIndex] == False) and (currentCost > costs[vertexIndex])):
                    currentCost = costs[vertexIndex]
                    currentVertex = vertexIndex
        self.printDijkstraResults(graph, costs, parents)


    # Void - FOR VALIDATION PURPOSES ONLY
    #        This function takes a MoGraph object and lists of the costs and parent vertices and edges
    #        constructed by running Dijkstra's algorithm on the MoGraph.  It steps through the parents
    #        list starting from the last mile marker and constructs a new reversed list of taken vertices.
    #        It saves the accumulated energy at the last mile marker vertex for comparison with
    #        the one-pass algorithm result, calls another function to convert (and to save) the reversed
    #        list of vertices to properly ordered stream tuples, and then prints out the tuples.
    #        NOTE:  Each element of takenVertices has a list of two elements, the start of the selected
    #               edge [0] and a reference to the edge [1].
    def printDijkstraResults(self, graph, costs, parents):
        lastVertexIndex = graph.nodeLookupDict[self.lastMileMarker]    # What index is the lastMileMarker vertex?
        totalCost = costs[lastVertexIndex]                             # Get the accumulated cost at the lastMileMarker
        takenVertices = [[self.lastMileMarker, None]]                  # Start building the taken vertex list from the end
        nextParent = lastVertexIndex                                   # Start looking at parent for the lastMileMarker
        firstVertexIndex = graph.nodeLookupDict[0]                     # Define loop to stop when we reach location 0
        completed = False                                              # - for loop completion test
        while((nextParent > -1) and (not completed)):                  # Loop while we have a defined parent and haven't reached 0
            if(parents[nextParent][0] > -1):                           # If valid parent index, append vertex and edge to list
                takenVertices.append(parents[nextParent][1:])
            if(nextParent == firstVertexIndex):                        # If we reached location 0, set loop completion
                completed = True
            else:                                                      # Otherwise, jump to next parent index
                nextParent = parents[nextParent][0]
        print 'Total minimum cost per Dijkstra is %i' % (totalCost)    # Print out cost and the number of vertices and edges
        print 'Graph has %i vertices and %i edges' % (graph.numVertices, graph.numEdges)
        self.dijkstraEnergy = totalCost                                # Save total cost for comparison to one-pass algorithm
        print self.convertVertexListToStreamTuples(takenVertices)      # Convert taken vertices to stream tuples, save, and print
        

    # Void - FOR VALIDATION PURPOSES ONLY
    #        This function takes a reversed list of taken vertices and edges and converts them to a 
    #        properly ordered stream tuple list.  It then saves the list for comparison with the
    #        one-pass algorithm result.
    #        NOTE:  Each element of takenVertices has a list of two elements, the start of the selected
    #               edge [0] and a reference to the edge [1].
    def convertVertexListToStreamTuples(self, takenVertices):
        streamTuples = []                               # Start with empty list at last vertex in the list
        lastVertex = len(takenVertices) - 1
        while(lastVertex > 0):                          # While last vertex > 0 (since we need to look at vertex - 1)
            start = takenVertices[lastVertex][0]        # Start is this vertex
            end = takenVertices[lastVertex - 1][0]      # End is the vertex preceding this one in the list
            startEdge = takenVertices[lastVertex][1]    # The edge object knows if it is a "default" edge or a "stream" edge
            if(do_debug):
                print 'DEBUG:  Converting vertex list: lastVertex %i, start %i, end %i' % (lastVertex, start, end)
                endEdge = takenVertices[lastVertex - 1][1]
                if(startEdge != None):
                    print 'DEBUG:  StartEdge start %i, end %i' % (startEdge.start, startEdge.end)
                else:
                    print 'DEBUG:  StartEdge -None-'
                if(endEdge != None):
                    print 'DEBUG:  EndEdge start %i, end %i' % (endEdge.start, endEdge.end)
                else:
                    print 'DEBUG:  EndEdge -None-'
            if(startEdge == None):                      # Error check - The start edge should never be None.
                print 'ERROR:  StartEdge is None at lastVertex %i' % (lastVertex)
            if(not startEdge.default):                  # If the edge isn't a "default" edge, then add it to the list of streams
                streamTuples.append((start, end))
            lastVertex = lastVertex - 1                 # Decrement to the next vertex and continue loop
        self.dijkstraTuples = streamTuples              # Save reference to list of tuples for comparison to one-pass algorithm
        return streamTuples                             # Return list reference to calling function


    # Void - FOR VALIDATION PURPOSES ONLY
    #        This function compares the one-pass and Dijkstra energy costs and stream tuple lists saved
    #        in the MoRev object.  It flags errors for all mismatches.
    #        If it finds a mismatch while stepping through the tuples, it will continue through the one-pass
    #        tuple results to see if any later tuples still match.
    def compareLinearToDijkstra(self):
        # Compare total energy costs and flag if inequal
        print 'Linear energy cost is %i, Dijkstra energy cost is %i' % (self.linearEnergy, self.dijkstraEnergy)
        if(self.linearEnergy != self.dijkstraEnergy):
            print 'ERROR:  The energy costs are not equal!'

        # Save and compare stream tuple list lengths; flag if inequal
        linearTuplesLength = len(self.linearTuples)
        dijkstraTuplesLength = len(self.dijkstraTuples)
        print 'Linear list length is %i, Dijkstra list length is %i' % (linearTuplesLength, dijkstraTuplesLength)
        if(linearTuplesLength != dijkstraTuplesLength):
            print 'ERROR:  The tuple lists should be the same length if they are equal!'

        # Compare all tuples in both lists
        allSame = True                              # start with no mismatch
        currentLinearTupleIndex = 0                 # start at one-pass tuple 0
        for currentDijkstraTupleIndex in range(0, dijkstraTuplesLength):  # over all Dijkstra tuples:
            if(currentLinearTupleIndex < linearTuplesLength):             # if current one-pass tuple index is still in range, check it
                currentLinearTuple = self.linearTuples[currentLinearTupleIndex]
                currentDijkstraTuple = self.dijkstraTuples[currentDijkstraTupleIndex]
                if((currentDijkstraTuple[0] == currentLinearTuple[0]) and   # check both start and end vertices; if match, print match
                   (currentDijkstraTuple[1] == currentLinearTuple[1])):
                    print 'cDT %i, cLT %i:  Both are (%i, %i)' % (currentDijkstraTupleIndex, currentLinearTupleIndex, currentDijkstraTuple[0],
                                                                  currentDijkstraTuple[1])
                    currentLinearTupleIndex = currentLinearTupleIndex + 1   # increment one-pass tuple index
                else:                                                       # otherwise, we found a mismatch
                    allSame = False                                         # set flag to indicate mismatch
                    savedLinearTupleIndex = currentLinearTupleIndex         # save the current one-pass tuple index then increment it
                    currentLinearTupleIndex = currentLinearTupleIndex + 1
                    foundMatch = False                                      # start loop with no match
                    while((not foundMatch) and (currentLinearTupleIndex < linearTuplesLength)):  # while no match and valid index
                        currentLinearTuple = self.linearTuples[currentLinearTupleIndex]          # check new one-pass tuple against Dijkstra
                        if((currentDijkstraTuple[0] == currentLinearTuple[0]) and       # if match, print match
                           (currentDijkstraTuple[1] == currentLinearTuple[1])):
                            print 'cDT %i, cLT %i:  Both are (%i, %i)' % (currentDijkstraTupleIndex, currentLinearTupleIndex,
                                                                          currentDijkstraTuple[0], currentDijkstraTuple[1])
                            foundMatch = True                                           # break loop by setting flag
                            currentLinearTupleIndex = currentLinearTupleIndex + 1       # increment current one-pass index for next loop
                    if(not foundMatch):                                     # otherwise, reset the one-pass index, print mismatch, loop
                        currentLinearTupleIndex = savedLinearTupleIndex
                        print 'cDT %i:  No match for tuple (%i, %i)' % (currentDijkstraTupleIndex, currentDijkstraTuple[0],
                                                                        currentDijkstraTuple[1])
            else:
                # We ran out of one-pass tuples to compare.  Print final statement and break out of loop.
                print 'Out of linear tuples.  cDT %i:  No match for tuple (%, %i)' % (currentDijkstraTupleIndex, currentDijkstraTuple[0],
                                                                                      currentDijkstraTuple[1])
                break
        if(allSame):   # Print success if all tuples matched
            print 'All tuples matched!'
        else:          # Else, print error
            print 'ERROR:  Not all tuples matched...'
        print

### END OF MoRev CLASS ###


################
##### GLOBAL FUNCTIONS
################

### Global function - (void) printUsage(doExit)
#
# Void - This function prints the expected script usage, including optional arguments for help, debug, and extra validation.
#        Pass the function a True value to exit the script after printing the usage.
def printUsage(doExit):
    if(len(sys.argv) > 0):
        scriptName = sys.argv[0]
    else:
        scriptName = 'Script'
    print scriptName + ' usage:  ' + scriptName + ' [-[d|h|v]] stream_input_file'
    print '  The script must be called with a stream input file.'
    print '  Optional arguments can be combined, can appear before or after the input file, and include:'
    print '    -d:  enable printing of extra debugging information; can be combined with -v'
    print '    -h:  print this help information and exit'
    print '    -v:  enable the very long, optional O(n^2) validation using an explicit graph'
    print
    if(doExit):
        exit(1)


### Global function - (String) parseArguments()
#
# String - This function parses the supplied arguments to the script and returns the file name in a string.
#          If an optional argument is incorrect or requests "help," then the function will print the script
#          usage and exit.  If more than one string is supplied without a leading dash, the function will
#          also print the script usage and exit.  It is flexible enough to support combined optional arguments.
def parseArguments():
    global do_debug, do_val_with_Dijkstra
    fileArgument = None              # Eventual holder of file name
    argLength = len(sys.argv)        # Number of arguments (0 is script name)
    if(argLength == 1):              # If only one item in the list, print usage and exit
        printUsage(True)
    for currentArgument in range(1, argLength):     # Over all arguments:
        argString = sys.argv[currentArgument]       # Get string and its length
        argStringLength = len(argString)            
        if(argStringLength > 0):                    # Check non-zero before indexing
            if(argString[0] == '-'):                # Is it an option and not a filename?
                if(argStringLength == 1):           # Just the dash?  Print usage and exit
                    printUsage(True)
                currentCharIndex = 1                # Otherwise, cycle through all options that may be in the string
                while(currentCharIndex < argStringLength):
                    if(argString[currentCharIndex] == 'h'):        # Is 'help' option? - Print usage and exit
                        printUsage(True)
                    if(argString[currentCharIndex] == 'v'):        # Is 'validate' option? - Turn on SLOW validation flag
                        do_val_with_Dijkstra = True
                    elif(argString[currentCharIndex] == 'd'):      # Is 'debug' option? - Turn on debug print statement flag
                        do_debug = True
                    else:
                        print '"' + argString[currentCharIndex] + '" is not a valid option.'
                        printUsage(True)
                    currentCharIndex = currentCharIndex + 1
            else:
                if(fileArgument != None):           # Do we already have a string that might be a file name?
                    print 'Please pass the script a maximum of one input file name.'
                    printUsage(True)                # If so, exit
                fileArgument = argString            # If not, save file name
    return fileArgument


################
##### MAIN 
################

fileArgument = parseArguments()       # Parse arguments, set any optional flags, and return file name
if(fileArgument == None):             # Exit if no file name has been supplied
    printUsage(True)
pathCalculator = MoRev(fileArgument)  # Instantiate the MoRev class object with the file; it runs the rest automatically



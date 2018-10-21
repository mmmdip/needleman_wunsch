import numpy as np
import sys

def read_file( filename ):
    '''
        reads the sequence file to get the sequences
    parameter:
        filename (str): name of the input file
    returns:
        seq (tuple): two separate strings/sequences 
    '''
    file = open( filename, 'r' )
    seq1 = file.readline().strip()
    seq2 = file.readline().strip()
    seq = ( seq1, seq2 )
    return seq

def build_dp_table( seq, param ):
    '''
        builds up the DP table using the parameters/scoring scheme passed
    parameters:
        seq (tuple): two strings/sequences
        param (tuple): the scoring scheme with the points/values
    returns:
        dp (list of lists): DP table with string/sequence alignment scores
        direction (list of lists): table with directions for alignment scores
    '''
    ( identity, transition, transversion, gap ) = param
    ( seq1, seq2 ) = seq
    seqLength = ( len( seq1 ) + 1, len( seq2 ) + 1 )
    
    # direction matrix is initiated with 'u's and 'l's representing up and left arrows in
    # the first row and first column
    direction = np.empty( seqLength, dtype = '<U10' )
    direction[ 1:, 0 ] = 'u'
    direction[ 0, 1: ] = 'l'
    
    dp = np.zeros( seqLength, dtype = int )
    
    # initialization of the dp_table
    for row in range( 1, seqLength[0] ):
        dp[row, 0] = dp[row - 1, 0] + gap
    for col in range( 1, seqLength[1] ):
        dp[0, col] = dp[0, col - 1] + gap
    
    for row in range( 1, seqLength[0] ):
        for col in range( 1, seqLength[1] ):
            # adding gap penalties for the values coming from left or up
            l = dp[row , col - 1] + gap
            u = dp[row - 1, col] + gap
            d = dp[row - 1, col - 1]
            char1 = seq1[row - 1]
            char2 = seq2[col - 1]
            if char1 == char2:
                # checking for match
                d += identity
            elif ( char1 in [ 'A', 'G' ] and char2 in [ 'A', 'G' ] ) or ( char1 in [ 'C', 'T' ] and char2 in [ 'C', 'T' ] ):
                # checking for transition
                d += transition
            else:
                # case for tranversion
                d += transversion
            score = max( l, u, d )
            dp[row, col] = score
            dir = ''
            if score == l:
                dir += 'l'
            if score == u:
                dir += 'u'
            if score == d:
                dir += 'd'
            # putting one or more arrow in the direction table:
            # l in string means left arrow
            # u in string means up arrow
            # d in string means diagonal arrow
            direction[row, col] = dir
    return ( dp, direction )

def trace_back( dp, direction, seq ):
    '''
        traces back the DP table to get the alignment/s
    parameters:
        dp (list of lists): DP table with string/sequence alignment scores
        direction (list of lists): table with directions for alignment scores
        seq (tuple): two strings/sequences
    returns:
        paths (list of lists): all possible sequence alignment paths
    '''
    ( seq1, seq2 ) = seq
    seqLength = ( len( seq1 ), len( seq2 ) )
    
    pos = ( seqLength[0], seqLength[1] )
    stk = []
    stk.append( pos )
    
    pathId = 0
    paths = [[]]
    
    divFlag = False
    while len( stk ) != 0:
        ( row, col ) = stk.pop()
        paths[pathId].append( ( row, col ))
        
        if direction[row, col] == 'u':
            stk.append( ( row - 1, col ) )
        elif direction[row, col] == 'd':
            stk.append( ( row - 1, col - 1 ) )
        elif direction[row, col] == 'l':
            stk.append( ( row, col - 1) )
        elif direction[row, col] == 'ld':
            stk.append( ( row - 1, col - 1 ) )
            stk.append( ( row, col - 1 ) )
            dividerNode = ( row, col )
        elif direction[row, col] == 'lu':
            stk.append( ( row - 1, col - 1 ) )
            stk.append( ( row, col - 1 ) )
            dividerNode = ( row, col )
            divFlag = True
        elif direction[row, col] == 'ud':
            stk.append( ( row - 1, col ) )
            stk.append( ( row - 1, col - 1 ) )
            dividerNode = ( row, col )
            divFlag = True
        elif direction[row, col] == 'lud':
            stk.append( ( row, col - 1 ) )
            stk.append( ( row - 1, col - 1 ) )
            stk.append( ( row - 1, col ) )
            dividerNode = ( row, col )
            divFlag = True
        elif direction[row, col] == '':
            paths.append([])
            pathId += 1
            if divFlag: paths[pathId] = paths[pathId - 1][:paths[pathId - 1].index( dividerNode ) + 1]
            continue
    paths = paths[:-1]
    return paths

def align_strings( paths, direction, seq ):
    '''
        aligns sequences based on the paths and directions passed as parameter
    paramters:
        paths (list of lists): all possible sequence alignment paths
        direction (list of lists): table with directions for alignment scores
        seq (tuple): two strings/sequences
    returns:
        alignedStrings (list of tuples): list containing the pair of aligned sequences
    '''
    alignedStrings = []
    ( seq1, seq2 ) = seq
    for path in paths:
        str1 = str2 = ''
        for id in range( 1, len( path )):
            ( curRow, curCol ) = path[id - 1]
            ( nextRow, nextCol ) = path[id]
            if nextRow + 1 == curRow and nextCol + 1 == curCol:
                # for diagonal arrow, both the characters from the sequences got added
                str1 += seq1[curRow - 1]
                str2 += seq2[curCol - 1]
            elif nextRow + 1 == curRow:
                # for left arrow, only the character from sequence 1 got added
                str1 += seq1[curRow - 1]
                str2 += '_'
            elif nextCol + 1 == curCol:
                # for up arrow, only the character from sequence 2 got added
                str1 += '_'
                str2 += seq2[curCol - 1]
        alignedStrings.append( ( str1[::-1], str2[::-1] ) )
    return alignedStrings
        
def get_best_alignment( alignedStrings, param ):
    '''
        calculates the best alignment from all the alignments based on the parameters passed
    parameters:
        alignedStrings (list of tuples): list containing the pair of aligned sequences
        param (tuple): the scoring scheme with the points/values
    returns:
        bestAlignment (tuple): the best sequence alignment calculated
        maxScore (int): the best alignment score
    '''
    scores = []
    ( identity, transition, transversion, gap ) = param
    
    for alignment in alignedStrings:
        ( align1, align2 ) = alignment
        score = 0
        for idx in range( len( align1 ) ):
            char1 = align1[idx]
            char2 = align2[idx]
            if char1 == char2:
                # for match, identity point is added to score
                score += identity
            elif ( char1 in [ 'A', 'G' ] and char2 in [ 'A', 'G' ] ) or ( char1 in [ 'C', 'T' ] and char2 in [ 'C', 'T' ] ):
                # for transition, transition penalty is deducted
                score += transition
            elif char1 == '_' or char2 == '_':
                # for gap, gap penalty is deducted
                score += gap
            else:
                # for transversion, transversion penalty is deducted
                score += transversion
        scores.append( score )
    ( bestAlignment, maxScore ) = ( alignedStrings[ scores.index( max( scores )) ], max( scores ) )
    return ( bestAlignment, maxScore )
            
def main():
    # getting the input filename from commandline
    filename = sys.argv[1]
    seq = read_file( filename )
    print( "1st input sequence:", seq[0] )
    print( "2nd input sequence:", seq[1] )
    
    identity = 4
    transition = -2
    transversion = -3
    gap = -8
    param = ( identity, transition, transversion, gap )
    
    # building dp_table
    ( dp, direction ) = build_dp_table( seq, param )
    
    # tracing back the dp_table for the alignments
    paths = trace_back( dp, direction, seq )
    
    # aligned strings are calculated for final ranking
    alignedStrings = align_strings( paths, direction, seq )
    
    # best alignment along with the best value is calculated
    ( bestAlignment, bestScore ) = get_best_alignment( alignedStrings, param )
    print( 'The best alignment is' )
    print( bestAlignment[0] )
    print( bestAlignment[1] )
    print( 'Best alignment score:', bestScore )
    
    dp_table = [ list( seq[0] ) ]
    for row in dp:
        dp_table.append( list( row ) )
    
    noOfArguments = len(sys.argv) - 1
    arguments = len(sys.argv) - 1

    position = 1  
    while (arguments >= position):  
        position = position + 1

    if sys.argv[noOfArguments] == '1':
        print( dp )  
    
if __name__ == '__main__':
    main()
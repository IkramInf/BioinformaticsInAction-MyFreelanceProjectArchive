import sys

# Dynamic programming solution to find optimal anchor sequence
def find_optimal_anchors(anchors):
    """
    Finds the optimal sequence of non-overlapping anchors using dynamic programming.

    Args:
        anchors (list): A list of anchor tuples, where each tuple contains two intervals.

    Returns:
        tuple: A tuple containing:
            - List of optimal anchors.
            - Maximum coverage achieved by the optimal anchors.
    """
    n = len(anchors)
    
    # Sort anchors by the ending position of the A and B intervals
    anchors.sort(key=lambda x: (x[0][1], x[1][1]))

    # DP table to store maximum span for each anchor
    dp = [0] * n
    # Table to reconstruct the optimal sequence
    prev = [-1] * n

    for i in range(n):
        # Calculate span of the current anchor
        span = anchors[i][0][1] - anchors[i][0][0]
        dp[i] = span  # Initialize DP value with the span

        # Find the last non-overlapping anchor
        prev_anchor = -1
        for j in range(i - 1, -1, -1):
            # Check for non-overlapping condition
            if anchors[j][0][1] < anchors[i][0][0] and anchors[j][1][1] < anchors[i][1][0]:
                prev_anchor = j
                break

        if prev_anchor != -1:
            # Update DP value with the maximum span including the previous anchor
            dp[i] = max(dp[i], dp[prev_anchor] + span)
            prev[i] = prev_anchor  # Update the previous anchor index

    # Find the index of the maximum coverage
    max_index = dp.index(max(dp))

    # Reconstruct the optimal sequence of anchors
    optimal_anchors = []
    while max_index != -1:
        optimal_anchors.append(anchors[max_index])
        max_index = prev[max_index]  # Move to the previous anchor

    optimal_anchors.reverse()  # Reverse to get the correct order
    return optimal_anchors, max(dp)  # Return the optimal anchors and their coverage

if __name__ == "__main__":
    # Check for the correct number of command-line arguments
    if len(sys.argv) != 2:
        print("Usage: python p7.py <alignment anchors>")
        sys.exit()

    # Read the alignment anchors from command-line argument
    alignment_anchors = sys.argv[1]  # Example: '0,5,0,5;6,10,7,11;4,8,4,8;12,25,16,29'
    
    anchors = []
    # Parse the alignment anchors into a list of tuples
    for anchor in alignment_anchors.split(';'):
        sA, eA, sB, eB = list(map(int, anchor.split(',')))  # Parse each anchor
        anchors.append(((sA, eA), (sB, eB)))  # Create a tuple of intervals

    # Find optimal anchors and their maximum coverage
    optimal_anchors, max_coverage = find_optimal_anchors(anchors)

    # Output the optimal anchors and the total coverage
    for optimal_anchor in optimal_anchors:
        print(f"{optimal_anchor[0]} -> {optimal_anchor[1]}")  # Print the optimal anchor intervals

    print("Coverage:", max_coverage)  # Print the total coverage

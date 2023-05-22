import sys

if len(sys.argv) != 5:
    print("usage: python recall.py <results_file> <target_col> <query_col> <num_true_positives>")
    exit(1)

results_file = sys.argv[1]
target_col = int(sys.argv[2])
query_col = int(sys.argv[3])
num_true_positives = int(sys.argv[4])

num_found = 0
with open(results_file) as file:
    for line in file:
        line_tokens = line.split()

        # the target name is formatted like:
        #     <target_name>/<id>/<from>-<to>, or
        #     <decoy#>
        target = line_tokens[target_col].split("/")[0]
        query = line_tokens[query_col]

        if target == query:
            num_found += 1
        elif line.startswith("decoy"):
            break

recall = num_found / num_true_positives

print(recall)
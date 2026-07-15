#!/bin/bash


pathing="./solvingProblems"

read -r -p "Enter lastLine file name: " lastLine

inputFile="$pathing/$lastLine"
sourceFile="$pathing/restartCreator.cpp"
executable="$pathing/restartCreator.exe"

# Make sure the CSV exists.
if [[ ! -f "$inputFile" ]]; then
    echo "Error: CSV file not found:"
    echo "$inputFile"
    exit 1
fi

# Compile restartCreator.cpp.
g++ "$sourceFile" -std=c++17 -o "$executable"

if [[ $? -ne 0 ]]; then
    echo "Error: restartCreator.cpp failed to compile."
    exit 1
fi

echo "Watching: $inputFile"
echo "Press Ctrl+C to stop manually."

# Number of seconds without a new completed line before stopping.
idle_limit=10
idle_seconds=0

count_complete_lines() {
    local count=0

    # This counts only lines ending with a newline.
    # A line that is still being written will not be counted yet.
    while IFS= read -r line; do
        ((count++))
    done < "$inputFile"

    printf '%s\n' "$count"
}

previous_count=$(count_complete_lines)

# Process the current final completed line once when starting.
if (( previous_count > 0 )); then
    echo "Processing current last line..."
    "$executable" "$inputFile"

    if [[ $? -ne 0 ]]; then
        echo "Warning: restartCreator.exe returned an error."
    fi
fi

while (( idle_seconds < idle_limit )); do
    sleep 1

    current_count=$(count_complete_lines)

    if (( current_count > previous_count )); then
        new_lines=$((current_count - previous_count))

        echo "$new_lines new completed line(s) detected."

        # Run restartCreator once for every newly completed line.
        for ((i = 0; i < new_lines; i++)); do
            "$executable" "$inputFile"

            if [[ $? -ne 0 ]]; then
                echo "Warning: restartCreator.exe returned an error."
            fi
        done

        previous_count=$current_count
        idle_seconds=0
    else
        ((idle_seconds++))
    fi
done

echo "The file has not received a completed line for $idle_limit seconds."
echo "Stopping file watcher."

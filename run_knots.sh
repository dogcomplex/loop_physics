#!/bin/bash
# Wrapper script to run knots.py and capture output

echo "Running knots.py and redirecting output to OUTPUT.txt..."

# Execute knots.py using the python interpreter in the current environment
# Redirect stdout (>) to OUTPUT.txt
# Redirect stderr (2>) to the same place as stdout (&1)
python knots.py > OUTPUT.txt 2>&1

# Check the exit code of the python script
EXIT_CODE=$?

if [ $EXIT_CODE -eq 0 ]; then
    echo "Finished running knots.py successfully. Check OUTPUT.txt for results."
else
    echo "Error running knots.py (Exit Code: $EXIT_CODE). Check OUTPUT.txt for details."
fi

exit $EXIT_CODE 
# Library to control keyboard and mouse ("pip install pyautogui" needed)
import pyautogui # type: ignore

# Library to measure time
import time

# Library to get the date and time
from datetime import datetime

# Library to get the folder path
import os

# Function to click in the position passed as parameters
def click_position(x, y):
    pyautogui.click(x, y)

# Function to autoclick in two places at the same time
def main():

    # Clicking Pylon Viewer
    click_position(224,67)
    t_start = time.perf_counter()

    # Clicking Waveforms
    click_position(1344,108)
    t_end = time.perf_counter()

    time_diff = (t_end-t_start)*1000
    
    # Time between each click
    print("Time between each line of code: " + str(time_diff) + " ms.")
    
    # Date and time to use as name of the file with the time between 
    # each click
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    name = f"time_diff_save_{timestamp}.txt"
    
    # Obtaining the folder path where we run this code
    current_root = os.path.dirname(os.path.realpath(__file__))
    
    # Name of the file + directory
    filename = os.path.join(current_root, name)
    
    # Saving the file in the same folder this code is
    with open(filename, "w") as file:
        file.write(str(time_diff))

# Call to main method
if __name__ == "__main__":
    main()
    print("Program finished.")
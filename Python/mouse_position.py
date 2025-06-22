# IMPORTANT - Sometimes the position printed is out of the bounds of the 
# individuals' monitor, in that case just restart the program or divide
# by 10, 100, 1000, ... to match something that is truthful.

# Library to control keyboard and mouse ("pip install pyautogui" needed)
import pyautogui

print("Move the mouse and look the command window to see coordinates " +
      "(press Ctrl+C to cancel).")

try:
    while True:
        
        # Position where is clicked is loaded to x and y variables
        x, y = pyautogui.position()
        
        # Position is printed and refreshed while moving the mouse
        print(f"Position: X={x}, Y={y}", end="\r")
        
except KeyboardInterrupt:
    print("\nProgram terminated.")
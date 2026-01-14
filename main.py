import tkinter as tk
from gui import RENDO

def main():
    try:
        root = tk.Tk()
        root.title("RENDO")
        app = RENDO(root)
        root.mainloop()
    except Exception as e:
        print("An error occurred while running the application:", str(e))

if __name__ == "__main__":
    main()


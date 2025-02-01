import pandas as pd
import matplotlib.pyplot as plt
import sys

def plot_loglog_order_of_convergence(file_path):
    """
    Generate a log-log plot of the order of convergence given a CSV file of errors.
    The plot will also show reference lines for first- and second-order convergence.
    
    Parameters:
    file_path (str): Path to the CSV file containing the error values.
    """
    # Load the CSV file
    data = pd.read_csv(file_path)
    
    # Extract the 'error' column as an array
    errors = data['Uerr'].values
    steps = data["Nodes"].values
    
    # Plot the error on a log-log scale
    plt.figure(figsize=(10, 6))
    plt.loglog(steps, errors, marker='o', linestyle='-', color='b', label='Error')
    
    # Add reference lines for first- and second-order convergence
    # Choose an initial reference point for the lines
    ref_point_x = steps[0]
    ref_point_y = errors[0]
    
    # First-order convergence reference line (slope = -1)
    first_order_y = ref_point_y * (steps / ref_point_x) ** -1
    plt.loglog(steps, first_order_y, 'r--', label='First-Order Convergence')
    
    # Second-order convergence reference line (slope = -2)
    second_order_y = ref_point_y * (steps / ref_point_x) ** -2
    plt.loglog(steps, second_order_y, 'g--', label='Second-Order Convergence')
    
    # Labels and title
    plt.xlabel('Step index (log scale)')
    plt.ylabel('Error (log scale)')
    plt.title('Log-Log Plot of Order of Convergence')
    plt.grid(True, which="both", ls="--")
    plt.legend()
    plt.savefig("conv.png")

# Example usage
# Ensure that 'output.csv' has a single column named 'error' with error values
if len(sys.argv) < 2:
    print("specify file path")
    exit()

csv_file_path = sys.argv[1]
plot_loglog_order_of_convergence(csv_file_path)

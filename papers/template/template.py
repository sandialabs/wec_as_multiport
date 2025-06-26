from fpdf import FPDF
import os

def create_pdf(filename, content):
    pdf = FPDF(format=(50, 50), unit='mm', orientation='L')
    pdf.add_page()
    pdf.set_font("Arial", size=48)
    pdf.rect(0, 0, pdf.w, pdf.h)
    pdf.multi_cell(0, 0, content, align='C')
    pdf.output(filename)

if __name__ == "__main__":
    
    out_dir = 'gfx'
    os.makedirs(out_dir, exist_ok=True)
    
    create_pdf(os.path.join(out_dir, "dummy.pdf"), "XX")
    
    # Get the current directory
    current_directory = os.getcwd()

    # List the contents of the current directory
    contents = os.listdir(current_directory)

    # Print the contents
    print(f"Current director: {current_directory}")
    print("Contents of the current directory:")
    for item in contents:
        print(item)

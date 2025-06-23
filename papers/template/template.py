# generate_pdfs.py
from fpdf import FPDF
import os

def create_pdf(filename, content):
    pdf = FPDF()
    pdf.add_page()
    pdf.rect(10, 10, pdf.w - 20, pdf.h - 20)
    pdf.set_font("Arial", size=24)
    pdf.cell(200, 10, txt=content, ln=True, align='C')
    pdf.output(filename)

if __name__ == "__main__":
    
    out_dir = 'gfx'
    os.makedirs(out_dir, exist_ok=True)
    
    create_pdf(os.path.join(out_dir, "template_figure1.pdf"), "This is Figure 1")
    create_pdf(os.path.join(out_dir, "template_figure2.pdf"), "This is Figure 2")
    
    # Get the current directory
    current_directory = os.getcwd()

    # List the contents of the current directory
    contents = os.listdir(current_directory)

    # Print the contents
    print(f"Current director: {current_directory}")
    print("Contents of the current directory:")
    for item in contents:
        print(item)


import sys

def cleanbib(lines_remove = ["abstract =", "urldate =", "file =", "note =", "url ="],       
            input_file_path = "bibliography.bib",
            output_file_path = "clean_bibliography.bib"):

    with open(input_file_path, "r", encoding="utf-8") as infile:
        lines = infile.readlines()

    cleaned_lines = [line for line in lines if not any(keyword in line for keyword in lines_remove)]

    with open(output_file_path, "w", encoding="utf-8") as outfile:
        outfile.writelines(cleaned_lines)

    print(f"Cleaned bibliography in {output_file_path}")

if __name__ == "__main__":
    if len(sys.argv) == 0:
        cleanbib()
    elif len(sys.argv) != 3:
        print("Usage: python3 cleanbib.py lines_remove bibliography.bib clean_bibliography.bib")
    else:
        cleanbib(sys.argv[1], sys.argv[2], sys.argv[3])


"""
usage: python3 cleanbib.py lines_remove bibliography.bib cleaned_bibliography.bib
"""
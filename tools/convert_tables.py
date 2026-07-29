import sys
import os
import argparse

def convert_file(input_path, output_path, var_name):
    with open(input_path, 'rb') as f:
        data = f.read()
    
    out_dir = os.path.dirname(output_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(output_path, 'w') as f_out:
        f_out.write("#pragma once\n")
        f_out.write("#include <array>\n\n")
        f_out.write("namespace CloudJ::Tables {\n")
        f_out.write(f"constexpr std::array<unsigned char, {len(data)}> {var_name} = {{\n")
        hex_data = ", ".join(f"0x{b:02x}" for b in data)
        f_out.write(hex_data)
        f_out.write("\n};\n")
        f_out.write("} // namespace CloudJ::Tables\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--test", action="store_true")
    parser.add_argument("--input-dir")
    parser.add_argument("--output-dir")
    args = parser.parse_args()

    if args.test:
        with open("test.dat", "wb") as f: f.write(b"test")
        convert_file("test.dat", "test.hpp", "test_var")
        os.remove("test.dat")
        os.remove("test.hpp")
        print("Verified")
        sys.exit(0)

    if args.input_dir and args.output_dir:
        for file in os.listdir(args.input_dir):
            if file.endswith(".dat"):
                var_name = file.replace(".dat", "").replace("-", "_")
                convert_file(os.path.join(args.input_dir, file), 
                             os.path.join(args.output_dir, f"{var_name}.hpp"), 
                             var_name)

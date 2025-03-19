# Specify Parent Directory
ROOT_DIR="Matrix"

# Specify Output Paths
VarX_outpath="$ROOT_DIR/VarX_Matrix.csv"
VarY_outpath="$ROOT_DIR/VarY_Matrix.csv"

# Process VarX
{ sed -e '1s/-a/Rep.1/g; 1s/-b/Rep.2/g; 1s/-c/Rep.3/g' "$ROOT_DIR/VarX_Header.csv"; } > "$VarX_outpath"
cat "$ROOT_DIR/VarX/"*.csv >> "$VarX_outpath"
# Sort the file
awk 'NR==1; NR > 1 {print | "sort -u -t, -k1,1"}' "$VarX_outpath" > "$VarX_outpath.tmp" && mv "$VarX_outpath.tmp" "$VarX_outpath"

# Process VarY
{ sed -e '1s/-a/Rep.1/g; 1s/-b/Rep.2/g; 1s/-c/Rep.3/g' "$ROOT_DIR/VarY_Header.csv"; } > "$VarY_outpath"
cat "$ROOT_DIR/VarY/"*.csv >> "$VarY_outpath"
# Sort the file
awk 'NR==1; NR > 1 {print | "sort -u -t, -k1,1"}' "$VarY_outpath" > "$VarY_outpath.tmp" && mv "$VarY_outpath.tmp" "$VarY_outpath"


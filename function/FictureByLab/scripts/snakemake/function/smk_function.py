import pandas as pd
import io
from Bio.bgzf import BgzfWriter

def get_coord_range(file):
    with open(file, 'r') as f:
        lines = f.readlines()
        xmin = float(lines[0].split()[1])
        xmax = float(lines[1].split()[1])
        ymin = float(lines[2].split()[1])
        ymax = float(lines[3].split()[1])
    
    offsetx = xmin
    rangex = int((xmax - xmin + 0.5) / 1 + 1)
    offsety = ymin
    rangey = int((ymax - ymin + 0.5) / 1 + 1)

    return offsetx, rangex, offsety, rangey

def sort_pixel_fun(coor,pixel_decode,bsize,scale,K,topk,output):
    # Get the coordinate range (offset and range) from the input coordinate
    offsetx, rangex, offsety, rangey = get_coord_range(coor)

    # Read the pixel decode file into a DataFrame
    df = pd.read_csv(pixel_decode, sep = '\t', skiprows = 1, header = None)
    df.columns = ['BLOCK', 'X', 'Y', 'K1', 'K2', 'K3', 'P1', 'P2', 'P3']


    df['BLOCK'] = ((df['X'] - offsetx) / bsize).astype(int) * bsize
    df['X'] = (((df['X'] - offsetx) * scale).astype(int)).clip(lower = 0)
    df['Y'] = (((df['Y'] - offsety) * scale).astype(int)).clip(lower = 0)


    df = df.sort_values(by = ['BLOCK', 'Y'])
    
    with BgzfWriter(f'{output}', 'wb') as f_out:
        f_out.write(f"##K={K};TOPK={topk}\n")
        f_out.write(f"##BLOCK_SIZE={bsize};BLOCK_AXIS=X;INDEX_AXIS=Y\n")
        f_out.write(f"##OFFSET_X={offsetx};OFFSET_Y={offsety};SIZE_X={rangex};SIZE_Y={rangey};SCALE={scale}\n")
        f_out.write(f"#BLOCK\tX\tY\tK1\tK2\tK3\tP1\tP2\tP3\n")
        

        for chunk in [df[i:i+10000] for i in range(0, df.shape[0], 10000)]:
            buffer = io.StringIO()
            chunk.to_csv(buffer, sep = '\t', index = False, header = False)
            buffer_bytes = buffer.getvalue().encode('utf-8')
            f_out.write(buffer_bytes)
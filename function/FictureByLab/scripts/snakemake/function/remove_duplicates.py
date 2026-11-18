import pandas as pd
import sys
import gzip
import io
from Bio.bgzf import BgzfWriter


def remove_duplicates(input, output):

    with gzip.open(input, 'rt') as f:  
        lines = f.readlines()
    
    data_start = next(i for i, line in enumerate(lines) if not line.startswith('#'))
    comments = [line for line in lines[:data_start]]

    df = pd.read_csv(input, sep='\t', compression='gzip', comment='#', header=data_start)
    df.columns = ['BLOCK', 'X', 'Y', 'K1', 'K2', 'K3', 'P1', 'P2', 'P3']

    df['original_order'] = df.index
    df_sorted = df.sort_values(by=['X', 'Y', 'P1', 'P2', 'P3'], ascending=[True, True, False, False, False])
    df_unique = df_sorted.drop_duplicates(subset=['X', 'Y'], keep='first')
    df_unique = df_unique.sort_values(by='original_order')
    df_unique = df_unique.drop(columns=['original_order'])
    
     
    with BgzfWriter(f'{output}', 'wb') as f_out:
        f_out.write(''.join(comments).encode('utf-8'))  

        for chunk in [df_unique[i:i+10000] for i in range(0, df_unique.shape[0], 10000)]:
            buffer = io.StringIO()
            chunk.to_csv(buffer, sep='\t', index=False, header=False)
            buffer_bytes = buffer.getvalue().encode('utf-8')
            f_out.write(buffer_bytes)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(1)
    input = sys.argv[1]
    output = sys.argv[2]
    remove_duplicates(input, output)

 
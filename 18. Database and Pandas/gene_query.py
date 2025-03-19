# Import the library
import argparse
import pandas as pd
import sqlalchemy

# section (a): Using the argparse module to add a ‘query’ argument
# Create the parser
parser = argparse.ArgumentParser()
# Add an argument
parser.add_argument('--query', type=str, required=True)
# Parse the argument
args = parser.parse_args()

# create engine for the provided sqlite database
arabidopsis = sqlalchemy.create_engine('sqlite:///arabidopsis.sqlite')

# create dataframe for the sqlite database
gene_annotation = pd.read_sql("select * from geneannotation", arabidopsis)

# keyword_search function to search query
def keyword_search(keyword):
    df = gene_annotation[gene_annotation.apply(lambda x: x.astype(str).str.contains(keyword).any(), axis=1)]
    return df

# section (b): Pass the query argument to the keyword_search function
match = keyword_search(args.query)

# section (c): Print the results dataframe
for index, gene in match.iterrows():
    print(gene['genename'], gene['annotation'])
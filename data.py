import pandas as pd

def load_data():

    df = pd.read_csv('F:\Proj\Gas_Analysis_Proj\data\dataset.csv')
    df = df.iloc[:, 1:]

    print(df)
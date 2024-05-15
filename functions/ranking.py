import pandas as pd
def filter_pregnancy_risk(data):
    dataPR = data[data["Pregnancy Risk"] != "No"]
    print("Data base with molecules dangerous for pregnant women")
    return dataPR
    

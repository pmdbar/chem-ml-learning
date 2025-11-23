import pandas as pd

def main():
    df = pd.DataFrame({
        "molecule": ["benzene", "toluene", "phenol"],
        "logP": [2.13, 2.73, 1.46]
    })

    print("Average logP:", df["logP"].mean())

if __name__ == "__main__":
    main()


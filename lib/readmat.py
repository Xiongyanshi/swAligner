# parse substitution matrix into nested dict

def readmat(filename):
    mat = {}
    ll = [i.split() for i in open(filename).readlines()]
    columns = ll[0][1:]
    columns = list(map(str.upper, columns))
    for row in ll[1:]:
        index = row[0].upper()
        mat[index] = {}
        for score, col in zip(row[1:], columns):
            mat[index][col] = float(score)
    return mat

def main():
    dnadefault = 'dna_default.mat'
    print(readmat(dnadefault))


if __name__ == '__main__':
    main()

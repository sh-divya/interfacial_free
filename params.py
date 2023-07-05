def ternary(c1, c2):
    return (0.5 * (c1 ** 0.5) + 0.5 * (c2 ** 0.5)) ** 2


if __name__ == "__main__":
    cmin_121 = 0.44
    cmin_313 = 0.759
    cmin_131 = 0.759
    cmin_232 = 0.49
    cmin_212 = 0.23
    cmin_313 = 0.759
    cmax_121 = 1.44
    cmax_313 = 2.8
    cmax_131 = 2.8
    cmax_232 = 2.8
    cmax_212 = 1.44
    cmax_313 = 2.8

    cmin_123 = ternary(cmin_121, cmin_313)
    cmax_123 = ternary(cmax_121, cmax_313)

    cmin_132 = ternary(cmin_131, cmin_232)
    cmax_132 = ternary(cmax_131, cmax_232)

    cmin_213 = ternary(cmin_212, cmin_313)
    cmax_213 = ternary(cmax_212, cmax_313)

    print("cmin_123, ", cmin_123)
    print("cmin_132, ", cmin_132)
    print("cmin_213, ", cmin_213)
    print("cmax_123, ", cmax_123)
    print("cmax_132, ", cmax_132)
    print("cmax_213, ", cmax_213)


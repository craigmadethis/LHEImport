import math
def convert_px_py_pz(px,py,pz):
    # transverse momentum
    pt = math.sqrt(math.fsum([math.pow(px, 2), math.pow(py, 2)]))
    # total momentum
    p = math.sqrt(math.fsum([math.pow(pt, 2), math.pow(pz, 2)]))
    if pt != 0:
        eta = math.asinh(pz / pt)
        phi = math.asin(py / pt)
    else:
        eta = math.copysign(float('inf'), pz)
        phi = 0
    return pt, eta, phi

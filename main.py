from hashlib import blake2b
from curves import Pallas
from fields import Fp

#     "iso-pallas",
A = Fp(0x18354a2eb0ea8c9c49be2d7258370742b74134581a27a59f92bb4b0b657a014b)
B = Fp(1265)
Z = Fp(-13)
theta = Fp(0x0f7bdb65814179b44647aef782d5cdc851f64fc4dc888857ca330bcc09ac318e)


def hash_to_field(curve_id: bytes, domain_prefix: bytes, message: bytes):
    suffix = domain_prefix + b'-' + curve_id + b'_XMD:BLAKE2b_SSWU_RO_' + \
             bytes([22 + len(curve_id) + len(domain_prefix)])
    hasher0 = blake2b(digest_size=64, person=b'\x00' * 16)
    hasher0.update(b'\x00' * 128 + message + b'\x00\x80\x00' + suffix)
    hasher1 = blake2b(digest_size=64, person=b'\x00' * 16)
    hasher1.update(hasher0.digest() + b'\x01' + suffix)
    hasher2 = blake2b(digest_size=64, person=b'\x00' * 16)
    hasher2.update(bytes(a ^ b for (a, b) in zip(hasher0.digest(), hasher1.digest())))
    hasher2.update(b'\x02' + suffix)
    element0 = Fp(int.from_bytes(hasher1.digest(), byteorder='big'))
    element1 = Fp(int.from_bytes(hasher2.digest(), byteorder='big'))
    return element0, element1


def map_to_curve_simple_swu(u):
    tv1 = Fp(1) / (Z ** 2 * u ** 4 + Z * u ** 2)  # 1. tv1 = inv0(Z^2 * u^4 + Z * u^2)
    x1 = ((Fp(0) - B) / A) * (Fp(1) + tv1)  # 2.  x1 = (-B / A) * (1 + tv1)
    if tv1 == Fp(0): x1 = B / (Z * A)  # 3.  If tv1 == 0, set x1 = B / (Z * A)
    gx1 = x1 ** 3 + A * x1 + B  # 4. gx1 = x1^3 + A * x1 + B
    x2 = Z * u ** 2 * x1  # 5.  x2 = Z * u^2 * x1
    gx2 = x2 ** 3 + A * x2 + B  # 6. gx2 = x2^3 + A * x2 + B
    if gx1.is_square():  # 7.  If is_square(gx1), set x = x1 and y = sqrt(gx1)
        x = x1
        y = gx1.sqrt()
    else:  # 8.  Else set x = x2 and y = sqrt(gx2)
        x = x2
        y = gx2.sqrt()
    if u.sgn0 != y.sgn0: y = Fp(0) - y  # 9.  If sgn0(u) != sgn0(y), set y = -y
    return x, y  # 10. return (x, y)


def is_on_iso_pallas_curve(pt):
    a = Fp(0x18354a2eb0ea8c9c49be2d7258370742b74134581a27a59f92bb4b0b657a014b)
    b = Fp(1265)
    return pt.z * pt.y ** 2 == pt.x ** 3 + a * pt.x * pt.z ** 2 + b * pt.z ** 3


def iso_map(pt):
    iso = [Fp(0x0e38e38e38e38e38e38e38e38e38e38e4081775473d8375b775f6034aaaaaaab),
           Fp(0x3509afd51872d88e267c7ffa51cf412a0f93b82ee4b994958cf863b02814fb76),
           Fp(0x17329b9ec525375398c7d7ac3d98fd13380af066cfeb6d690eb64faef37ea4f7),
           Fp(0x1c71c71c71c71c71c71c71c71c71c71c8102eea8e7b06eb6eebec06955555580),
           Fp(0x1d572e7ddc099cff5a607fcce0494a799c434ac1c96b6980c47f2ab668bcd71f),
           Fp(0x325669becaecd5d11d13bf2a7f22b105b4abf9fb9a1fc81c2aa3af1eae5b6604),
           Fp(0x1a12f684bda12f684bda12f684bda12f7642b01ad461bad25ad985b5e38e38e4),
           Fp(0x1a84d7ea8c396c47133e3ffd28e7a09507c9dc17725cca4ac67c31d8140a7dbb),
           Fp(0x3fb98ff0d2ddcadd303216cce1db9ff11765e924f745937802e2be87d225b234),
           Fp(0x025ed097b425ed097b425ed097b425ed0ac03e8e134eb3e493e53ab371c71c4f),
           Fp(0x0c02c5bcca0e6b7f0790bfb3506defb65941a3a4a97aa1b35a28279b1d1b42ae),
           Fp(0x17033d3c60c68173573b3d7f7d681310d976bbfabbc5661d4d90ab820b12320a),
           Fp(0x40000000000000000000000000000000224698fc094cf91b992d30ecfffffde5)]

    x = pt.x / pt.z
    y = pt.y / pt.z
    x_num = iso[0] * x ** 3 + iso[1] * x ** 2 + iso[2] * x + iso[3]
    x_den = x ** 2 + iso[4] * x + iso[5]
    y_num = iso[6] * x ** 3 + iso[7] * x ** 2 + iso[8] * x + iso[9]
    y_den = x ** 3 + iso[10] * x ** 2 + iso[11] * x + iso[12]
    return type(pt)(x_num / x_den, (y * y_num) / y_den, Fp(1))


if __name__ == '__main__':
    tmp1, tmp2 = hash_to_field(b'pallas', b'z.cash:test', b'Trans rights now!')
    (act_x1, act_y1) = map_to_curve_simple_swu(tmp1)
    assert act_x1 == Fp(0x05c3482fe40155e152fdc0be06c4766b67a2b3d8d9bb64ee6137382879dc2160)
    assert act_y1 == Fp(0x3825fb730c259375175ff31b94dc36dcf031b13f3116bda725f1c98717739f1f)
    (act_x2, act_y2) = map_to_curve_simple_swu(tmp2)
    assert act_x2 == Fp(0x2c6e5aa1a88cd76c8a9d436438d2993244bf7704e4f322a86d0890bd6cee28ab)
    assert act_y2 == Fp(0x0b20c46efea44d15e4828808c86a72789d54328635ba4274d8e9b48d9654f65b)

    q0 = Pallas(act_x1, act_y1, Fp(1))
    q1 = Pallas(act_x2, act_y2, Fp(1))
    r = q0 + q1
    assert r == Pallas(Fp(0x3da8497a87f06e28b7983f044f5f93575daf4806e0735700ebd79184070bb58e),
                       Fp(0x271b2b52a5e759ee28a21db1a520739b1f53a1960433d08593ec10e225dec8c0),
                       Fp(1))
    print("r is: ", r)

    assert is_on_iso_pallas_curve(r)

    z = iso_map(r)
    assert z == Pallas(Fp(0x1818cda31ffdc8c3ff23df3d88c26f952340257d0f187a0236695c9b640b6bd3),
                       Fp(0x1e20888510123752166a0306332e126289f6f9a2774160395f2f1efc9b1280c),
                       Fp(1))
    print("iso out: ", z)

import des
import operator
from functools import reduce
import sympy
import random
import os
import math

def _sub(i, S_BOX):
    # replicate SBOX behavior of DES
    # block is an array of 6-bits
    # we return an array of 4-bits
    block = bin(i)[2:].zfill(6)
    ret = []
    row = int(str(block[0])+str(block[5]), 2)#Get the row with the first and last bit
    column = int(''.join([str(x) for x in block[1:][:-1]]),2) #Column is the 2,3,4,5th bits
    val = S_BOX[row][column] #Take the value in the SBOX appropriated for the round (i)
    for j in range(3,-1,-1):
        ret.append((val >> j) & 1)            
    return int(''.join(map(str, ret)), 2) 

def xor_bits(x):
    # xor all bits of a integer x
    return reduce(operator.xor, map(int, bin(x)[2:]))

def generate_LAT(S_BOX):
    Lat = [[-32 for i in range(15)] for j in range(63)]
    for alpha in range(1, 64): # input mask from 1 to 63
        for beta in range(1, 16): # output mask from 1 to 16
            for x in range(64):
                if xor_bits(x & alpha) == xor_bits(_sub(x, S_BOX) & beta): # check if the input bits xorred equals the output bits xorred after substitution and increase the bias in the LAT if it does
                    Lat[alpha - 1][beta - 1] += 1
    return Lat

def find_linear_approx(lat_idx, THRESHOLD=14, ignore_indexes=None):
    # search lat and print indxes with abs(bias) above THRESHOLD
    # by default THRESHOLD is 14
    # if ignore_indexes is set, it's should be a list containing values from 0 to 7 and if those bit indexes are used in the mask then we can't use them
    # important for when charset is not truely "random" and certain bits are biased, we can only use bits that are 100% biased or 0% biased.
    lat = LATs[lat_idx]
    if ignore_indexes:
         for i in range(len(lat)):
            for j in range(len(lat[0])):
                if abs(lat[i][j]) >= THRESHOLD:
                    input_mask = i + 1
                    output_mask = j + 1
                    a = [des.E[lat_idx*6 + x] for x, y in enumerate(bin(input_mask)[2:].zfill(6)) if y == "1"] # input bits
                    # output bits after sbox and xor is permuted using the P array
                    b = [des.P.index(lat_idx*4 + x) for x, y in enumerate(bin(output_mask)[2:].zfill(4)) if y == "1"] # output bits
                    if any(thing % 8 in ignore_indexes for thing in a + b): continue
                    print(f"input mask: {bin(i + 1)[2:].zfill(6)}, output mask: {bin(j + 1)[2:].zfill(4)}, bias: {lat[i][j]}")
    else:
        for i in range(len(lat)):
            for j in range(len(lat[0])):
                if abs(lat[i][j]) >= THRESHOLD:
                    print(f"input mask: {bin(i + 1)[2:].zfill(6)}, output mask: {bin(j + 1)[2:].zfill(4)}, bias: {lat[i][j]}")

LATs = []
for i in range(8):
    LATs.append(generate_LAT(des.S_BOX[i]))

#for i in range(8):
#    print(f"SBOX #{i + 1}:")
#    find_linear_approx(i, THRESHOLD=8, ignore_indexes=[3, 4, 6])
#    print("--------------------------")

'''
SBOX #1:
input mask: 010000, output mask: 1111, bias: -18
input mask: 111011, output mask: 1000, bias: -14
input mask: 111111, output mask: 0001, bias: -14
--------------------------
SBOX #2:
input mask: 011111, output mask: 1000, bias: -14
input mask: 100010, output mask: 1011, bias: -16
input mask: 100110, output mask: 0100, bias: -14
input mask: 110111, output mask: 0100, bias: -14
--------------------------
SBOX #3:
input mask: 010010, output mask: 1101, bias: -14
input mask: 100010, output mask: 1111, bias: 16
input mask: 101111, output mask: 1000, bias: 14
input mask: 111011, output mask: 0001, bias: -14
--------------------------
SBOX #4:
input mask: 100010, output mask: 1111, bias: -16
input mask: 101000, output mask: 1111, bias: -16
input mask: 101011, output mask: 0110, bias: 16
input mask: 101011, output mask: 1001, bias: -16
--------------------------
SBOX #5:
input mask: 010000, output mask: 1111, bias: -20
input mask: 100010, output mask: 1110, bias: -16
input mask: 100100, output mask: 1101, bias: 14
input mask: 111111, output mask: 0100, bias: 14
--------------------------
SBOX #6:
input mask: 010000, output mask: 0111, bias: -14
input mask: 100010, output mask: 1011, bias: 14
--------------------------
SBOX #7:
input mask: 010000, output mask: 1111, bias: -14
input mask: 010010, output mask: 1011, bias: -14
input mask: 100010, output mask: 1110, bias: -16
input mask: 101000, output mask: 1101, bias: 14
input mask: 111011, output mask: 0100, bias: -18
input mask: 111111, output mask: 0001, bias: 14
--------------------------
SBOX #8:
input mask: 000100, output mask: 1111, bias: 14
input mask: 010000, output mask: 1111, bias: -16
input mask: 100010, output mask: 1110, bias: 16
input mask: 101100, output mask: 0111, bias: -14
--------------------------
'''

# Picking this linear approx for SBOX 5
# input mask: 010000, output mask: 1111, bias: -20

'''
R_expand = [R[x] for x in E]
tmp = xor(self.round_keys[i], R_expand)
tmp = self._substitute(tmp)
# we linear approx _substitute
# for example we pick S_5 with an input mask of 16 and an output mask of 15
# S1 => 0 to 3, S2 => 4 to 7 (bit indexes)
# S5 => 16 to 19
# then these bits are permuted by P
# output mask of 15 means we use all bits from 16 to 19
P = [14, 13, 16, 26, 2, 19, 5, 22, 30, 27, 24, 0, 11, 9, 29, 3, 8, 15, 25, 12, 7, 23, 17, 21, 6, 28, 18, 4, 20, 10, 31, 1]
so bits 16 to 19 in x are now shuffled to indexes 2, 5, 22, and 26, which is our new  indexes for the output of round function F
tmp = [tmp[x] for x in P]
L, R = R, xor(L, tmp)
'''

# This function prints out the bits related to the equation for the linear approx we pick for a specific SBOX
def print_matsui_equation(sbox_num, input_mask, output_mask):
    # We represent the round function F of DES using our linear approximation
    # 32-bit input is expanded to 48 using the E array so we have to apply this transformation per what sbox and mask we're using
    a = [des.E[(sbox_num - 1)*6 + x] for x, y in enumerate(bin(input_mask)[2:].zfill(6)) if y == "1"] # input bits
    # output bits after sbox and xor is permuted using the P array
    b = [des.P.index((sbox_num - 1)*4 + x) for x, y in enumerate(bin(output_mask)[2:].zfill(4)) if y == "1"] # output bits
    c = [(sbox_num - 1)*6 + x for x, y in enumerate(bin(input_mask)[2:].zfill(6)) if y == "1"] # key bits
    a, b, c = sorted(a), sorted(b), sorted(c)

    print(f"X[{','.join(map(str, a))}] ^ F(X, K)[{','.join(map(str, b))}] = K[{','.join(map(str, c))}]")

#print_matsui_equation(5, 0b010000, 0b1111)

# SBOX 5
# input mask: 001000, output mask: 1011, bias: -12
# print_matsui_equation(5, 0b001000, 0b1011)
# X[17] ^ F(X, K)[2,5,26] = K[26]

# SBOX 6
# input mask: 001000, output mask: 0111, bias: -8
#print_matsui_equation(6, 0b001000, 0b0111)
# X[21] ^ F(X, K)[7,21,23] = K[32]

# SBOX 3
# input mask: 010000, output mask: 0110, bias: 8
#print_matsui_equation(3, 0b010000, 0b0110)
# X[8] ^ F(X, K)[13,29] = K[13]

# X[16] ^ F(X, K)[2,5,22,26] = K[25]
# above says we're using bit 25 of the key xorred with bit 16 of our input
# and it's strongly linearly related to bits 2, 5, 22, and 26 of the output of the round function F

# Due to the Feistal construction of DES, we have to pick a second linear approximation to cancel out this one
# It has to use only a singular bit index 16 for the round function F
# For now assume it's only a single input index
def find_linear_approx_chain(inp_bit_idx, THRESHOLD=10, ignore_indexes=None):
    # (sbox_num - 1)*4 + x == P[inp_bit_idx]
    # 0 <= sbox_num <= 7
    # 0 <= x <= 3
    inp_bit_idx = des.P[inp_bit_idx]
    sbox_num = inp_bit_idx // 4
    o = inp_bit_idx % 4
    output_mask = '0'*o + '1' + '0'*(3 - o) # assuming single index
    print(f"SBOX #{sbox_num + 1}:")
    print(f"Output mask: {output_mask}")

    output_mask = int(output_mask, 2)
    lat = LATs[sbox_num]
    # print possible input masks
    # similar to before we're just gonna print out the ones that are abs(bias) >= THRESHOLD
    # we lower the default THRESHOLD here since there's a lot less choices
    for i in range(len(lat)):
        if abs(lat[i][output_mask - 1]) >= THRESHOLD:
            if ignore_indexes:
                input_mask = i + 1
                a = [des.E[sbox_num*6 + x] for x, y in enumerate(bin(input_mask)[2:].zfill(6)) if y == "1"] # input bits
                # output bits after sbox and xor is permuted using the P array
                b = [des.P.index(sbox_num*4 + x) for x, y in enumerate(bin(output_mask)[2:].zfill(4)) if y == "1"] # output bits
                if any(thing % 8 in ignore_indexes for thing in a + b): continue
                print(f"input mask: {bin(i + 1)[2:].zfill(6)}, output mask: {bin(output_mask)[2:].zfill(4)}, bias: {lat[i][output_mask - 1]}")
            else:
                print(f"input mask: {bin(i + 1)[2:].zfill(6)}, output mask: {bin(output_mask)[2:].zfill(4)}, bias: {lat[i][output_mask - 1]}")

# find_linear_approx_chain(16)
#find_linear_approx_chain(17, THRESHOLD=0, ignore_indexes=[3, 4, 6])
#find_linear_approx_chain(21, THRESHOLD=0, ignore_indexes=[3, 4, 6])
#find_linear_approx_chain(8, THRESHOLD=0, ignore_indexes=[3, 4, 6])

""" 
SBOX #3:
Output mask: 1000
input mask: 101111, output mask: 1000, bias: 14
input mask: 110111, output mask: 1000, bias: 10
input mask: 111111, output mask: 1000, bias: -12 
"""

# We can pick input mask: 101111, output mask: 1000, bias: 14
# Now we have 2 chains

# print_matsui_equation(5, 0b010000, 0b1111)
# X[16] ^ F(X, K)[2,5,22,26] = K[25]

# print_matsui_equation(3, 0b101111, 0b1000)
# X[7,9,10,11,12] ^ F(X, K)[16] = K[12,14,15,16,17]

# Using matsui's construction we use the first approx for rounds 2/4 and the second one for rounds 1/5

def xor_multiple(arr, bits):
    out = 0
    for bit in bits:
        out = arr[bit] ^ out
    return out

K = sympy.IndexedBase('K') # key
PH = sympy.IndexedBase('PH') # upper half of input
PL = sympy.IndexedBase('PL') # lower half of input
CH = sympy.IndexedBase('CH') # upper half of output
CL = sympy.IndexedBase('CL') # lower half of output
X = sympy.IndexedBase('X') # intermediate input to found function F

# represent equation as left and right side
# if an expression already exists, then remove it, since xor of the same thing is 0

def k(rnd_num, key_bits):
    return [K[rnd_num, key_bit] for key_bit in key_bits]

def round_func(rnd_num, input_bits, output_bits, key_bits):
    # X[input_bits] ^ F(X, K)[output_bits] = K[key_bits]
    left = [*r(rnd_num, input_bits), *l(rnd_num, output_bits), *r(rnd_num+1, output_bits)]
    right = k(rnd_num, key_bits)
    return left, right

def l(rnd_num, bit_idxes):
    if rnd_num == 1: return [PH[bit_idx] for bit_idx in bit_idxes]
    return r(rnd_num - 1, bit_idxes)

def r(rnd_num, bit_idxes):
    if rnd_num == 1: return [PL[bit_idx] for bit_idx in bit_idxes]
    if rnd_num == 5: return [CL[bit_idx] for bit_idx in bit_idxes]
    if rnd_num == 6: return [CH[bit_idx] for bit_idx in bit_idxes]
    return x(rnd_num, bit_idxes)

def x(rnd_num, bit_idxes):
    return [X[rnd_num, bit_idx] for bit_idx in bit_idxes] 


def remove_duplicates(l, r):
    # remove elements that show up even amount of times
    l_new = []
    for element in l:
        if l.count(element) % 2 != 0:
            l_new.append(element)
    
    r_new = []
    for element in r:
        if r.count(element) % 2 != 0:
            r_new.append(element)

    return l_new, r_new


# AB-BA
#l1, r1 = round_func(1, [7,9,10,11,12], [16], [12,14,15,16,17])
#l2, r2 = round_func(2, [16], [2,5,22,26], [25])
#l4, r4 = round_func(4, [16], [2,5,22,26], [25])
#l5, r5 = round_func(5, [7,9,10,11,12], [16], [12,14,15,16,17])
#left, right = remove_duplicates([*l1, *l2, *l4, *l5], [*r1, *r2, *r4, *r5])
#print(left)
#print(right)

'''
[PL[7], PL[9], PL[10], PL[11], PL[12], PH[16], PL[2], PL[5], PL[22], PL[26], CL[2], CL[5], CL[22], CL[26], CL[7], CL[9], CL[10], CL[11], CL[12], CH[16]]
[K[1, 12], K[1, 14], K[1, 15], K[1, 16], K[1, 17], K[2, 25], K[4, 25], K[5, 12], K[5, 14], K[5, 15], K[5, 16], K[5, 17]]
'''
# Above equation matches up with matsui's analysis

# X[17] ^ F(X, K)[2,5,26] = K[26]
# print_matsui_equation(4, 0b001000, 0b0001)
# X[13] ^ F(X, K)[17] = K[20]

l1, r1 = round_func(1, [13], [17], [20])
l2, r2 = round_func(2, [17], [2,5,26], [26])
l4, r4 = round_func(4, [17], [2,5,26], [26])
l5, r5 = round_func(5, [13], [17], [20])
left, right = remove_duplicates([*l1, *l2, *l4, *l5], [*r1, *r2, *r4, *r5])
print(left)
print(right)

'''
[PL[13], PH[17], PL[2], PL[5], PL[26], CL[2], CL[5], CL[26], CL[13], CH[17]]
[K[1, 20], K[2, 26], K[4, 26], K[5, 20]]
'''


def algorithm_1(pairs, bit_idxs, bias):
    # given an equation
    # P[i1, i2, ..., in] ^ C[j1, j2, ..., jb] = K[k1, k2, ..., kc]
    # We have N plaintexts and calculate when left side of the equation is 0
    # Let's say that happens T times
    # if T > N / 2
    # then K[k1, k2, ..., kc] = 0 (when p > 0.5) otherwise 1
    # else
    # K[k1, k2, ..., kc] = 1 (when p > 0.5) otherwise 0
    # p is the probability that our linear approximation holds
    # for example a bias of -20 means the linear approx holds with probability 12/64

    N = len(pairs)
    ph_bits, pl_bits, ch_bits, cl_bits = bit_idxs

    cnt = 0
    for _pt, _ct in pairs:
        pt = des.bytes_to_bits(_pt)
        ct = des.bytes_to_bits(_ct)
        pt = [pt[x] for x in des.IP] # apply initial permutation
        ct = [ct[x] for x in des.inverse_FP] # apply inverse of final permutation

        ph, pl = pt[:32], pt[32:]
        ch, cl = ct[:32], ct[32:]

        T = xor_multiple(ph, ph_bits) ^ xor_multiple(pl, pl_bits) ^ xor_multiple(ch, ch_bits) ^ xor_multiple(cl, cl_bits)
        if T == 0:
            cnt += 1
    
    if cnt > N // 2:
        if bias > 0: return 0, cnt
        return 1, cnt
    elif bias > 0: return 1, cnt
    return 0, cnt

N = 65536

bias1 = (-12 + 32)/32
bias2 = (4 + 32)/32
prob = bias1*bias2 + (1 - bias1)*(1 - bias2)
total_prob = prob**2 + (1 - prob)**2
print("Probability Linear Approximation holds:", total_prob) # probability overall linear approx equation holds
print(f"# of plaintexts needed: {math.ceil(abs(total_prob - 0.5)**(-2))}")


# TEST ALGORITHM 1
for _ in range(64):
    key = os.urandom(8)
    flag = "WACON2023{" + key.hex() + "}"
    print(f"Trial #{_ + 1}")
    print(f"flag = {flag}")
    key_bits = des.bytes_to_bits(key)

    tmp_key = des.bits_to_bytes(key_bits)
    cipher = des.DES(tmp_key, 5)
    rkeys = cipher.round_keys

    k1, k2, k4, k5 = rkeys[0], rkeys[1], rkeys[3], rkeys[4]

    random.seed(os.urandom(32))
    pairs = []
    for _ in range(N):
        pt = ''.join(random.choices('0123456789-', k=8)).encode()
        #pt = os.urandom(8)
        ct = cipher.encrypt(pt)
        pairs.append((pt, ct))

    guess_key_bit, cnt = algorithm_1(pairs, [[17], [13,2,5,26], [17], [13,2,5,26]], total_prob)
    print(f"Count: {cnt}, Bias: {(cnt - N//2)/(N//2)}")
    print("Guessed key bit:", guess_key_bit)
    print("Correct key bit:", xor_multiple(k1, [20]) ^ xor_multiple(k5, [20]) ^ k2[26] ^ k4[26])
    print("------------------------")

""" charset = '0123456789-'
for i in range(8):
    cnt = 0
    for c in charset:
        b = (ord(c) >> i) & 1
        if b == 0: cnt += 1
    print(f"Bit #{7 - i}: Bias {cnt/len(charset) - 0.5}") """
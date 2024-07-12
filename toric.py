#!/usr/bin/env sage
#from sage.all import *
from PIL import Image, ImageDraw
import sys
import math
from ortools.linear_solver import pywraplp

def main():
    if '-s' in sys.argv:
        print("mosaic string?:")
        mosaic_string = input()
        toric_mosaic.string_catalog(mosaic_string) 
        return
    if '-r' in sys.argv: 
        print("p?:", end = " ")
        p = int(input())
        print("q?:", end = " ")
        q = int(input())
        
        if math.gcd(p,q) != 1:
            print("p and q must be relatively prime")
            return
        toric_mosaic.rapunzel_mosaic(p,q)
        return
    #Prints if no methods were called
    print("Usage: python3 toric.py [option]")
    print("Valid options:")
    print("\t-h: print this message")
    print("\t-r: Produce toric mosaic of torus knot (p,q) using rapunzel algorithm.")
    print("\t-s: Catalog mosaic from string. (not fully implemented)")

class toric_mosaic:
    valid_connections = (
        (()),
        ((2,3),(3,2)),
        ((0,3),(3,0)),
        ((0,1),(1,0)),
        ((2,1),(1,2)),
        ((2,0),(0,2)),
        ((1,3),(3,1)),
        ((2,3),(3,2),(1,0),(0,1)),
        ((0,3),(3,0),(2,1),(1,2)),
        ((0,2),(1,3),(2,0),(3,1)),
        ((0,2),(1,3),(2,0),(3,1))
    )
    @classmethod
    # Written by Luc Ta, https://github.com/luc-ta/torus-knot-toric-mosaics
    def get_rapunzel_params(cls, p,q):
        # Create the mip solver with the SCIP backend.
        solver = pywraplp.Solver.CreateSolver("SAT")
        if not solver:
            raise Exception("solver issue")

        infinity = solver.infinity()
        # Make h and v non-negative integer variables
        h = solver.IntVar(0.0, infinity, "h")
        v = solver.IntVar(0.0, infinity, "v")

        # First constraint
        solver.Add(-3 * h - v - p + q + 4 >= 0)

        # Second constraint: r >= -R
        solver.Add(q - 2 * (h + v + p) + 4 >= 0)
        
        # Third constraint: R >= r
        solver.Add(h + 3 * v <= q - 3 * p + 4)
        
        # Fourth constraint
        solver.Add(h >= v)

        # Optimize n, the size of the toric mosaic.
        solver.Minimize(q - h - v)

        status = solver.Solve()

        if status == pywraplp.Solver.OPTIMAL:
            return [h.solution_value(),v.solution_value()]
        else:
            raise Exception("No Optimal Parameters Found")


    @classmethod
    def rapunzel_mosaic(cls, p, q):
        mosaic = []
        try:
            rapunzel_params = cls.get_rapunzel_params(p,q)
        except Exception as exp:
            print(f"Error: {exp.args}")
            return
        h = int(rapunzel_params[0])
        v = int(rapunzel_params[1])
        r = (p+v) - h if v != 0 else 2 - h
        n = q - h - v
        print(f"n = {n}\nh = {h}\nv = {v}")
        for i in range(0,p-2):
            mosaic.extend([7]*n)
        for i in range(0,h):
            row = [7]*i + [10]*(p-1)
            row.extend([7]*(n-len(row)))
            mosaic.extend(row)
        if v != 0:
            for i in range(0,(p+v-2)):
                row = [8]*(h-1-i) + [9]*(min(i + 1,p-1,(p+v-2)-i))        
                row.extend([8]*(n-len(row)))
                mosaic.extend(row)
        if r < 0:
            for i in range(0,abs(r)):
                mosaic.extend([8]*n)
        if r > 0:
            for i in range(0, abs(r)):
                mosaic.extend([7]*n)
        while len(mosaic) < n**2:
            mosaic.extend([6]*(n))
        
        for i in range(0, len(mosaic)):
            if i % n == 0:
                print("")
            print(f"{mosaic[i]:x}", end = "")

        to_png(mosaic, f"images/{p},{q}:{n}n,{h}h,{v}v.png")


    @classmethod
    def string_catalog(cls, mosaic_string):
        #initializing things here instead of wasting time with allocation between batches
        size = int(len(mosaic_string)**(0.5))
        mosaic = [9]*(2*(size ** 2))
        satisfied = [False]*(2*(size ** 2))
        crossing_number = [0]*(2*(size ** 2))
        gauss_code = []
        crossing_signs = []
        made_connections = [[] for _ in range(size ** 2)]
        crossing_count = 0
        i = 0
        num = 0
        curr_tile = 0
        starting_tile = None
        face = 0

        for char in mosaic_string:
            num = int(char, base = 10)
            mosaic[i] = num
            satisfied[i] = num == 0
            if starting_tile == None and num != 0:
                starting_tile = i
            i += 1

        curr_tile = starting_tile
        face = cls.valid_connections[mosaic[curr_tile]][0][0] 
        while curr_tile != starting_tile or not satisfied[curr_tile]:
            for conn in cls.valid_connections[mosaic[curr_tile]]:
                if conn[0] == face:
                    made_connections[curr_tile].append(conn)
                    if ((len(made_connections[curr_tile]) == 1) and mosaic[curr_tile] < 7) or (len(made_connections[curr_tile]) == 2):
                        satisfied[curr_tile] = True
                    if conn in [(0,3),(1,2)]:
                        down_cusps += 1
                    if conn in [(3,0),(2,1)]:    
                        up_cusps += 1

                    #Dealing with crossings
                    if mosaic[curr_tile] == 9:
                        if satisfied[curr_tile]:
                            if conn[0] % 2 == 1:
                                gauss_code.append(crossing_number[curr_tile])
                            else:
                                gauss_code.append(-crossing_number[curr_tile])

                            if conn[0] + made_connections[curr_tile][0][0] == 3: #Positive crossing -- works because the sum of starting faces for a positive crossing is always 3
                                crossing_signs[crossing_number[curr_tile]-1] = 1
                            else:
                                crossing_signs[crossing_number[curr_tile]-1] = -1
                        else:
                            crossing_count += 1
                            crossing_signs.append(0)
                            if conn[0] % 2 == 1:
                                crossing_number[curr_tile] = crossing_count
                                gauss_code.append(crossing_count)
                            else:
                                crossing_number[curr_tile] = crossing_count
                                gauss_code.append(-crossing_count)

                    if mosaic[curr_tile] == 10:
                        if satisfied[curr_tile]:
                            if conn[0] % 2 == 0:
                                gauss_code.append(crossing_number[curr_tile])
                            else:
                                gauss_code.append(-crossing_number[curr_tile])

                            if conn[0] + made_connections[curr_tile][0][0] == 3: #Positive crossing -- works because the sum of starting faces for a positive crossing is always 3
                                crossing_signs[crossing_number[curr_tile]-1] = 1
                            else:
                                crossing_signs[crossing_number[curr_tile]-1] = -1
                        else:
                            crossing_count += 1
                            crossing_signs.append(0)
                            if conn[0] % 2 == 1:
                                crossing_number[curr_tile] = crossing_count
                                gauss_code.append(crossing_count)
                            else:
                                crossing_number[curr_tile] = crossing_count
                                gauss_code.append(-crossing_count)
                    face = (conn[1] + 2) % 4
                    if face == 0: #left
                        curr_tile -= 1
                    elif face == 1: #down
                        curr_tile += size
                    elif face == 2: #right
                        curr_tile += 1
                    elif face == 3: #up 
                        curr_tile -= size
                    break

        #if all(satisfied) and len(gauss_code) != 0: 
            #knot_type = Link([[gauss_code],crossing_signs]).get_knotinfo()
        #else:
            #knot_type = 'idk'
        #to_png(mosaic,f"{knot_type}|{writhe - (up_cusps + down_cusps // 2)}|{abs(up_cusps - down_cusps) // 2}.png")
        print(f"{gauss_code},{crossing_signs}")
            #tup = ('0_1' if len(gauss_code) == 0 else (cls.knot_list.get(Link([[gauss_code],crossing_signs]).homfly_polynomial(), Link([[gauss_code],crossing_signs]).homfly_polynomial())), writhe - (up_cusps + down_cusps) // 2, abs(up_cusps - down_cusps) // 2)
            #out_file.write(f"{'0_1' if len(gauss_code) == 0 else cls.knot_list[Link([[gauss_code],crossing_signs]).homfly_polynomial()]} | {writhe - (up_cusps + down_cusps) // 2} | {abs(up_cusps - down_cusps) // 2} | {mosaic_string}\n")

def to_png(matrix,output_filename):
    tile_size = 64
    border_size = 4
    border_color = (196, 196, 196, 255)
    size = int(len(matrix)**0.5)

    tile_images = {}
    for num in range(11):
        file_name = f"tiles/{num}.png"
        try:
            tile_images[num] = Image.open(file_name).convert("RGBA")
        except FileNotFoundError:
            print(f"Failed to load image {file_name}")

    mosaic_width = size * tile_size + 2 * border_size
    mosaic = Image.new("RGBA", (mosaic_width, mosaic_width), border_color)
    draw = ImageDraw.Draw(mosaic)

    for i, tile in enumerate(matrix):
            if tile in tile_images:
                img_tile = tile_images[tile]
                for y in range(tile_size):
                    for x in range(tile_size):
                        pixel = img_tile.getpixel((x, y))
                        mosaic.putpixel(( (i % size) * tile_size + x + border_size, (i // size) * tile_size + y + border_size), pixel)

    mosaic.save(output_filename)


main()

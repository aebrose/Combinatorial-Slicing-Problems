#Sage hyperplane arrangement gets error if you have redundant hyperplanes 
#only necessary for central slicing arrangment
def delete_redundant_normals(verts):
    verts_no_reps = []
    for vertex in verts:
        add_vertex = true
        for s in verts_no_reps:
            if Matrix([vector(QQ, vertex), vector(QQ, s)]).rank() == 1:
                add_vertex = false
                break
        if add_vertex == true:
            verts_no_reps.append(vertex)
    return verts_no_reps


def central_hyperplane_arrangement(P):
    normals = delete_redundant_normals(P.vertices())
    hyperplane_data = [[list(n), 0] for n in normals]
    H = HyperplaneArrangements(QQ, tuple(['x'+str(i) for i in range(1, P.ambient_dim()+1)])) 
    return H(hyperplane_data)

def refined_central_hyperplane_arrangement(P,Q):
    normals = delete_redundant_normals(P.vertices_list() + Q.vertices_list())
    hyperplane_data = [[list(n), 0] for n in normals]
    H = HyperplaneArrangements(QQ, tuple(['x'+str(i) for i in range(1, P.ambient_dim()+1)]))
    return H(hyperplane_data)

def check_pointwise_cse(P, Q, verbose=False):
    H_PQ = refined_central_hyperplane_arrangement(P,Q)
    points = sample_normalvectors(H_PQ, all_or_generic="all")
    cse = True
    for pt in points:
        alpha_u = [0]+list(pt)
        H = Polyhedron(eqns=[alpha_u])
        P_slice = P.intersection(H)
        Q_slice = Q.intersection(H)
        if not P_slice.is_combinatorially_isomorphic(Q_slice):
            cse = False
            if verbose:
                print(f"\nDIFF SLICES FOR  H=Poly(eqns = {alpha_u}).\n We print f-vectors of the slices (maybe enough).")
                print(f"\nP_slice f-vector: {P_slice.f_vector()}")
                print(f"\nQ_slice f-vector: {Q_slice.f_vector()}")
            return cse 
    return cse 

def sample_normalvectors(hp_arr, all_or_generic="all"):
    """
    Input: hyperplane arrangement.
    Output: one sample vector from each chamber/face.

    For "generic": samples one vector from each full-dimensional region.
    For "all": samples one vector from each face of each region, with duplicates removed.
    note: If we extend to unbounded polyhedra than do something similar but combine .rays() with .lines() too.
    """
    regions = hp_arr.regions()
    normalvectors = []        

    for reg in regions: #full-dimensional regions
        rays = reg.rays() #if fails then add if len(rays)==0: continue
        n = sum([vector(list(r)) for r in rays]) #try sum(vector(QQ, list(r)) for r in rays) 
        if -n not in normalvectors:
            normalvectors.append(n)
        
            if all_or_generic == "all": #include lower-dimensional regions 
                for k in range(1, reg.dimension()): # skip empty cone
                    for face_of_region in reg.faces(k):
                        rays = face_of_region.as_polyhedron().rays()
                        n = sum([vector(list(r)) for r in rays])
                        if n not in normalvectors and -n not in normalvectors: #cause we will do many faces double
                            normalvectors.append(n)
    return normalvectors


def combtype_is_already_in(list_polytopes, S):
    """
    Input: a list of polyhedra (as Sage objects), and another polytope S (think: slices of a polytope)
    Output: true or false, depending on whether combinatorial type of S is in list of polytopes or not
    """
    b = false
    i = 0
    while(b == false and i < len(list_polytopes)):
        b = S.is_combinatorially_isomorphic(list_polytopes[i])
        i += 1
    return b

def find_central_slices(P, all_or_generic = "all", set_or_multiset = "set"):
    """
    Input: polytope P, parameter options "all" or "generic, and "set" or "multiset"
    Output: List of slice polytopes (as Sage objects)
    """
    H_P = central_hyperplane_arrangement(P)
    points = sample_normalvectors(H_P, all_or_generic)
    central_slices = []
    
    for pt in points:
        H = Polyhedron(eqns=[[0]+list(pt)])
        Q = P.intersection(H)
    
        if set_or_multiset == "set" and (combtype_is_already_in(central_slices, Q) == false):
            central_slices.append(Q)
        if set_or_multiset == "multiset":
            central_slices.append(Q)
        
    return [S for S in central_slices if S.dim()==P.dim()-1]

def find_affine_slices(P, all_or_generic = "all", set_or_multiset = "set"):
    verts = P.vertices()
    new_verts = []
    
    for vertex in verts:
        new_vertex = vertex + [1]
        new_verts.append(new_vertex)
    
    P_new = Polyhedron(vertices = new_verts)
    slices_on_height1 = find_central_slices(P_new, all_or_generic, set_or_multiset)
    slices = []
    for S in slices_on_height1:
        verts_h1 = S.vertices_list()
        verts = [v[:-1] for v in verts_h1]
        slices.append(Polyhedron(vertices = verts))
    
    return [S for S in slices if S.dim() == P.dim()-1]

def plot_poly(P, poly_color="green", poly_opacity = 0.5, edge_color = "black", vertex_color="black", w_origin = True, plot_it = True, viewpoint = ([-0.5,-0.5,-0.7],105)):
    G = Graphics()
    G += P.plot(alpha=poly_opacity, color=poly_color)
    for e in P.faces(1):
        G += Polyhedron(e).plot(point=False, frame=False, alpha=0.5, color=edge_color, thickness=1.5)
    for v in P.faces(0):
        G += Polyhedron(v).plot(frame=False, alpha=0.9, color=vertex_color, size=11)
    if w_origin:
        G += Polyhedron(vertices = [[0,0,0]]).plot(frame=False, alpha=1, color="black", size=22)
    if plot_it:
        G.show(viewpoint = viewpoint, online=True)
    return G 

def plot_poly_w_slice(P, P_slice_reps, poly_color="teal", slice_color="blue", w_origin = False, viewpoint = ([-0.5,-0.5,-0.7],105)):
    list_of_graphics = []
    for S in P_slice_reps:
        if S.dim() == P.dim() - 1:
            print(f"{S.n_vertices()}-gon slice")
            G = (plot_poly(P, poly_color=poly_color,  poly_opacity=0.2, plot_it = False, w_origin=w_origin) + plot_poly(S, poly_color=slice_color, edge_color = slice_color, vertex_color=slice_color, poly_opacity=0.4, plot_it = False, w_origin=w_origin))
            G.show(viewpoint = viewpoint, online=True)
            list_of_graphics.append(G)
    return list_of_graphics

print("""
        ============================================================
        Algorithm to compute slices — FUNCTIONS
        ============================================================
        1. find_central_slices(P, all_or_generic="all", set_or_multiset="set")
                options: "all" or "generic", "set" or "multiset"
                Output: list of Sage polyhedra objects representing a slice per category
                    (all/set, all/multiset, generic/set, generic/multiset)

        2. find_affine_slices(P, all_or_generic="all", set_or_multiset="set")
                options: as above
                Output: list of Sage polyhedra objects representing a slice per category

        3. plot_poly(P, poly_color="green", poly_opacity = 0.5,
                edge_color = "black", vertex_color = "black", w_origin = True,
                plot_it = True, viewpoint = ([-0.5,-0.5,-0.7],105))
                    
        4. plot_poly_w_slice(P, P_slice_reps, poly_color="teal",
                slice_color="blue", w_origin = False, viewpoint = ([-0.5,-0.5,-0.7],105))

        5.check_pointwise_cse(P, Q, verbose=False)
                output: whether P cap u^perp ~ Q cap u^perp for all u 
        ============================================================
        """)

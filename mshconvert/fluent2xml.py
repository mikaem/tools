#!/usr/bin python
import re
from sets import Set
from copy import copy
from numpy import zeros, array
from dolfin import File, MeshEditor, Mesh, Cell, facets

def fluent2xml(ifilename, ofilename):
    """Converting from ANSYS Fluent format (.msh) to FEniCS xml format
    The fluent mesh (the .msh file) is basically stored as a list of vertices, and then a 
    list of faces for each zone of the mesh, the interior and the boundaries."""    
    
    # Use regular expressions to identify sections and tokens found in a fluent file
    re_dimline   = re.compile(r"\(2\s(\d)\)")
    re_comment   = re.compile(r"\(0\s.*")
    re_zone_init = re.compile(r"\(10\s\(0\s(\w+)\s(\w+)\s(\d+)\s(\d+)\)\)")
    re_zone      = re.compile(r"\(10\s\((\w+)\s(\w+)\s(\w+)\s(\d+)\s(\d)\)(\(|)")
    re_face_init = re.compile(r"\(13(\s*)\(0\s+(\w+)\s+(\w+)\s+(0|0 0)\)\)")
    re_face      = re.compile(r"\(13(\s*)\((\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\)(\s*)(\(|)")
    re_periodic  = re.compile(r"\(18.*\((\w+)\s+(\w+)\s+(\w+)\s+(\w+)\).*\(")
    re_pfaces    = re.compile(r"((^\s)|)(\w+)(\s*)(\w+)")
    re_cells_init= re.compile(r"\(12(\s*)\(0(\s+)(\w+)(\s+)(\w+)(\s+)(0|0 0)\)\)")
    re_cells     = re.compile(r"\(12.*\((\w+)\s+(\w+)\s+(\w+)\s+(\d+)\s+(\d+)\)\)")
    re_cells2    = re.compile(r"\(12(\s*)\((\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\)(\s*)(\(|)")
    re_zones     = re.compile(r"\((45|39)\s+\((\d+)\s+(\S+)\s+(\S+).*\)\((.*|[0-9]+[\.]*[0-9]*)\)\)")
    re_parthesis = re.compile(r"(^\s*\)(\s*)|^\s*\)\)(\s*)|^\s*\(\s*)")

    # Declare som maps that will be built when reading in the lists of vertices and faces:
    cell_map = {}               # Maps cell id with vertices
    boundary_cells = {}         # List of cells attached to a boundary facet. Key is zone id
    zones = {}                  # zone information (not really used yet)

    def read_periodic(ifile, periodic_dx):
        """Scan past periodic section. Periodicity is computed by FEniCS."""
        while 1:
            line = ifile.readline()
            a = re.search(re_pfaces, line)
            if a:
                continue
            break

    def read_zone_vertices(dim, Nmin, Nmax, ifile, editor):
        """Scan ifile for vertices and add to mesh_editor."""
        # First line could be either just "(" or a regular vertex.
        # Check for initial paranthesis. If paranthesis then read a new line, else reset
        pos = ifile.tell()
        line = ifile.readline()
        if not re.search(re_parthesis, line): 
            ifile.seek(pos) # reset            
        # read Nmax-Nmin vertices
        for i in range(Nmin, Nmax + 1):
            line = ifile.readline()
            vertex = [eval(x) for x in line.split()]
            if dim == 2:
                editor.add_vertex(i - Nmin, vertex[0], vertex[1])
            else:
                editor.add_vertex(i - Nmin, vertex[0], vertex[1], vertex[2])
        
    def read_faces(zone_id, Nmin, Nmax, bc_type, face, ifile):
        """Read all faces and create cell_map + boundary maps."""        
        pos = ifile.tell() # current position
        line = ifile.readline()
        if not re.search(re_parthesis, line): # check for initial paranthesis. If paranthesis then read a new line, else reset
            ifile.seek(pos)
            
        # read Nmax-Nmin faces
        for i in range(Nmin, Nmax + 1):
            line = ifile.readline()
            ln = line.split()
            if face == 0:
                nd = int(ln[0], 16) # Number of vertices
                nds = [int(x, 16) for x in ln[1:(nd + 1)]]
                cells = [int(x, 16) for x in ln[(nd + 1):]]
            else:
                nd = face
                nds = [int(x, 16) for x in ln[:nd]]
                cells = [int(x, 16) for x in ln[nd:]]
                        
            if min(cells) == 0: # A boundary zone
                if zone_id in boundary_cells:
                    boundary_cells[zone_id][max(cells)] = array(nds)
                else:
                    boundary_cells[zone_id] = {max(cells): array(nds)}

            for c in cells:
                if c > 0:    
                    if not c in cell_map:
                        cell_map[c] = copy(nds)                    
                    else:
                        cell_map[c] = list(Set(cell_map[c] + nds))
                                        
    def scan_fluent_mesh(ifile, mesh, editor):  
        """Scan fluent mesh and generate maps."""
        dim = 0
        one = 0
        while 1:
            line = ifile.readline()
            if len(line) == 0:
                print 'Finished reading file\n'
                break

            if dim == 0: # Dimension usually comes first
                a = re.search(re_dimline, line)
                if a: 
                    print 'Reading dimensions\n'
                    dim = int(a.group(1))
                    editor.open(mesh, dim, dim)
                    continue
            
            if one == 0: # The total number of vertices
                a = re.search(re_zone_init, line)
                if a:
                    print 'Reading zone info\n'
                    one, num_vertices, dummy1, dummy2 = int(a.group(1)), \
                        int(a.group(2), 16), int(a.group(3), 16), int(a.group(4))
                    editor.init_vertices(num_vertices)
                    continue
                
            a = re.search(re_zone, line) # Vertices
            if a:
                zone_id, first_id, last_id = int(a.group(1), 16), \
                    int(a.group(2), 16), int(a.group(3), 16)
                print 'Reading ', last_id - first_id + 1,' vertices in zone ', zone_id + 1, '\n'
                read_zone_vertices(dim, first_id, last_id, ifile, editor)
                continue
                
            a = re.search(re_zones,line) # Zone info
            if a:
                print 'Reading zone info ', line
                dummy, zone_id, zone_type, zone_name, radius =  \
                        int(a.group(1)), int(a.group(2)),  a.group(3), \
                        a.group(4), a.group(5)
                zones[zone_id] = [zone_type, zone_name, radius]
                continue
            
            a = re.search(re_cells_init, line) # Get total number of cells/elements
            if a:
                print 'Reading cell info ', line
                first_id, tot_num_cells = int(a.group(3),16), int(a.group(5), 16)
                editor.init_cells(tot_num_cells)
                continue

            a = re.search(re_cells,line) # Get the cell info.
            if a:
                zone_id, first_id, last_id, bc_type, element_type = \
                    int(a.group(1)), int(a.group(2), 16), int(a.group(3), 16), \
                    int(a.group(4), 16), int(a.group(5), 16)
                print 'Found ', last_id - first_id + 1,' cells in zone ', zone_id, '\n'
                if last_id == 0:
                    raise TypeError("Zero elements!")
                continue

            a = re.search(re_cells2,line) # Get the cell info.
            if a:
                raise TypeError("Wrong cell type. Can only handle one single cell type")

            a = re.search(re_face_init, line)
            if a:
                print 'Reading total number of faces\n', line
                continue
                
            a = re.search(re_face, line)
            if a:
                print 'Reading faces ', line
                zone_id, first_id, last_id, bc_type, face_type = \
                    int(a.group(2), 16), int(a.group(3), 16), int(a.group(4), 16), \
                    int(a.group(5), 16), int(a.group(6), 16)
                read_faces(zone_id, first_id, last_id, bc_type, face_type, ifile)
                continue
            
            a = re.search(re_periodic, line)
            if a:
                print 'Scanning past periodic connectivity\n', line
                read_periodic(ifile, periodic_dx)
                continue
            
            if any([re.search(st, line) for st in (re_parthesis, re_comment)]) or \
                                                                not line.strip():
                continue
                    
            # Should not make it here
            raise IOError('Something went wrong reading fluent mesh.')        

    def write_fenics_file(ofile, mesh, editor):
        
        dim = mesh.geometry().dim()
        for i in range(1, len(cell_map)+1):
            if dim == 2:
                editor.add_cell(i-1, cell_map[i][0]-1, cell_map[i][1]-1, cell_map[i][2]-1)
            else:
                editor.add_cell(i-1, cell_map[i][0]-1, cell_map[i][1]-1, cell_map[i][2]-1, cell_map[i][3]-1)
        
        mesh.order()
        # Set MeshValueCollections from info in  boundary_cell
        mvc = mesh.domains().markers(dim-1)
        for zone, cells in boundary_cells.iteritems():
            for cell, nds in cells.iteritems():
                dolfin_cell = Cell(mesh, cell-1)
                vertices_of_cell = dolfin_cell.entities(0)
                vertices_of_face = nds - 1
                for jj, ff in enumerate(facets(dolfin_cell)):
                    facet_vertices = ff.entities(0)
                    if all(map(lambda x: x in vertices_of_face, facet_vertices)):
                        local_index = jj
                        break
                mvc.set_value(cell-1, local_index, zone)
            
        ofile << mesh        
        from dolfin import plot
        plot(mesh, interactive=True)
        print 'Finished writing FEniCS mesh\n'

    ifile  = open(ifilename, "r")   
    ofile  = File(ofilename)
    mesh = Mesh()
    editor = MeshEditor()
    scan_fluent_mesh(ifile, mesh, editor)
    write_fenics_file(ofile, mesh, editor)
    ifile.close()

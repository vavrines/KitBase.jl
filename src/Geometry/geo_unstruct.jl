"""
    struct UnstructPSpace{A,B,C,D,E,F,G,H,I,J,K,L} <: AbstractPhysicalSpace
        cells::A # all information: cell, line, vertex
        points::B # locations of vertex points
        cellid::C # node indices of elements
        cellType::D # inner/boundary cell
        cellNeighbors::E # neighboring cells id
        cellFaces::F # cell edges id
        cellCenter::G # cell center location
        cellArea::H # cell size
        facePoints::I # ids of two points at edge
        faceCells::J # ids of two cells around edge
        faceCenter::K # edge center location
        faceType::L
    end

Physical space with unstructured mesh
"""
struct UnstructPSpace{A,B,C,D,E,F,G,H,I,J,K,L} <: AbstractPhysicalSpace
    cells::A # all information: cell, line, vertex
    points::B # locations of vertex points
    cellid::C # node indices of elements
    cellType::D # inner/boundary cell
    cellNeighbors::E # neighboring cells id
    cellFaces::F # cell edges id
    cellCenter::G # cell center location
    cellArea::H # cell size
    facePoints::I # ids of two points at edge
    faceCells::J # ids of two cells around edge
    faceCenter::K # edge center location
    faceType::L
end

function UnstructPSpace(file::T) where {T<:AbstractString}
    cells, points = read_mesh(file)
    cellid = extract_cell(cells)
    edgePoints, edgeCells, cellNeighbors = mesh_connectivity_2D(cellid)
    cellType = mesh_cell_type(cellNeighbors)
    cellArea = mesh_area_2D(points, cellid)
    cellCenter = mesh_center_2D(points, cellid)
    edgeCenter = mesh_face_center(points, edgePoints)
    cellEdges = mesh_cell_face(cellid, edgeCells)
    edgeType = mesh_face_type(edgeCells, cellType)

    return UnstructPSpace(
        cells,
        points,
        cellid,
        cellType,
        cellNeighbors,
        cellEdges,
        cellCenter,
        cellArea,
        edgePoints,
        edgeCells,
        edgeCenter,
        edgeType,
    )
end


function write_vtk(ks::T1, ctr) where {T1<:AbstractSolverSet}
    KitBase.write_vtk(ks.ps, ctr)

    return nothing
end

function KitBase.write_vtk(ps::UnstructPSpace, ctr)
    cdata = zeros(length(ctr), length(ctr[1].w))
    for i in eachindex(ctr)
        cdata[i, :] .= ctr[i].prim
        if size(cdata, 2) > 1
            cdata[i, end] = 1.0 / cdata[i, end]
        end
    end
    KitBase.write_vtk(ps.points, ps.cellid, cdata)

    return nothing
end

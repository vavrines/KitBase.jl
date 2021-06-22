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
struct UnstructPSpace{
    A,
    B<:AbstractMatrix{<:AbstractFloat},
    C<:AbstractMatrix{<:Integer},
    D<:AbstractVector{<:Integer},
    E<:AbstractVector{<:AbstractFloat},
    F<:AbstractArray{<:AbstractFloat,3},
} <: AbstractPhysicalSpace
    cells::A # all information: cell, line, vertex
    points::B # locations of vertex points
    cellid::C # node indices of elements
    cellType::D # inner/boundary cell
    cellNeighbors::C # neighboring cells id
    cellFaces::C # cell edges id
    cellCenter::B # cell center location
    cellArea::E # cell size
    cellNormals::F # cell unit normal vectors
    facePoints::C # ids of two points at edge
    faceCells::C # ids of two cells around edge
    faceCenter::B # edge center location
    faceType::D # inner/boundary face
    faceArea::E # face area
end

function UnstructPSpace(file::T) where {T<:AbstractString}
    cells, points = read_mesh(file)
    p = mesh_connectivity_2D(cells, points)

    return UnstructPSpace((cells, points)..., p...)
end


function write_vtk(ks::T1, ctr) where {T1<:AbstractSolverSet}
    write_vtk(ks.ps, ctr)

    return nothing
end

function write_vtk(ps::UnstructPSpace, ctr)
    cdata = zeros(length(ctr), length(ctr[1].w))
    for i in eachindex(ctr)
        cdata[i, :] .= ctr[i].prim
        if size(cdata, 2) > 1
            cdata[i, end] = 1.0 / cdata[i, end]
        end
    end
    write_vtk(ps.points, ps.cellid, cdata)

    return nothing
end

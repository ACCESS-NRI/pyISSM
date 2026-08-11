"""
Context: `gmshplanet(md, radius, resolution, refine=<prior mesh>,
refinemetric=<per-vertex array>)` needs to size the new mesh according to
values defined on the *old* mesh's vertices. The Implementation plan's
original sketch (`gmsh.view.addModelData(view_tag, 0, 'planet', 'NodeData',
node_tags, ...)`) assumed the prior mesh's node tags already existed in the
new 'planet' model - they don't, since 'planet' is a fresh OCC sphere with
no mesh yet at that point.

Resolution, prototyped and verified here: Gmsh Views are session-global, not
tied to whichever model is "current" when later used as a background field.
So the prior mesh's nodes/triangles are loaded as a *real* mesh into their
own throwaway discrete model, the View's NodeData is attached to that model,
and the 'PostView' background field on the new model references the View by
tag alone - no shared node tags required between old and new meshes.

Note: `Mesh.MeshSizeMin`/`Mesh.MeshSizeMax` are
session-global options. If they were set narrow for the initial (unrefined)
mesh pass and never widened, they silently clamp the background field to a
uniform size on the refined pass, making refinement appear to do nothing.
"""
import numpy as np


def attach_refine_background_field(gmsh, node_tags, node_coords,
                                    tri_tags, tri_node_tags, metric,
                                    source_model_name='refine_source'):
    """
    Load a prior mesh's triangulation + per-vertex metric into a throwaway
    discrete Gmsh model and return a field tag ready for
    `gmsh.model.mesh.field.setAsBackgroundMesh(field_tag)` on a *different*,
    currently-active model.

    Parameters
    ----------
    gmsh : module
        The initialized `gmsh` module (session must already be `gmsh.initialize()`d).
    node_tags : (nv,) int array
        Prior mesh's Gmsh node tags, as returned by `gmsh.model.mesh.getNodes()`.
    node_coords : (3*nv,) float array
        Prior mesh's flat node coordinates, as returned by `getNodes()`.
    tri_tags : (nt,) int array
        Prior mesh's triangle element tags (type-2 elements only).
    tri_node_tags : (3*nt,) int array
        Prior mesh's flat triangle node-tag connectivity.
    metric : (nv,) float array
        Desired local element size at each prior-mesh vertex, same order as `node_tags`.
    source_model_name : str
        Name for the throwaway model that hosts the transferred mesh; must
        not collide with an existing model name in the session.

    Returns
    -------
    field_tag : int
        A 'PostView' mesh-size field tag. Caller sets it as the active
        model's background mesh; does not switch the active model itself.
    """
    metric = np.asarray(metric, dtype=float)
    if metric.shape[0] != len(node_tags):
        raise ValueError(
            f'metric length {metric.shape[0]} != node_tags length {len(node_tags)}')

    calling_model = gmsh.model.getCurrent()

    gmsh.model.add(source_model_name)
    gmsh.model.setCurrent(source_model_name)
    surf_tag = gmsh.model.addDiscreteEntity(2)
    gmsh.model.mesh.addNodes(2, surf_tag, node_tags, node_coords)
    gmsh.model.mesh.addElements(2, surf_tag, [2], [tri_tags], [tri_node_tags])

    view_tag = gmsh.view.add('refinemetric')
    gmsh.view.addModelData(
        view_tag, 0, source_model_name, 'NodeData',
        node_tags, metric.reshape(-1, 1))

    gmsh.model.setCurrent(calling_model)

    field_tag = gmsh.model.mesh.field.add('PostView')
    gmsh.model.mesh.field.setNumber(field_tag, 'ViewTag', view_tag)
    return field_tag


if __name__ == '__main__':
    # Standalone demonstration, run directly: `python scripts/gmsh_refine_transfer.py`
    import gmsh

    RADIUS = 1000.0

    gmsh.initialize()
    gmsh.option.setNumber('General.Terminal', 0)

    gmsh.model.add('prior')
    gmsh.model.occ.addSphere(0, 0, 0, RADIUS)
    gmsh.model.occ.synchronize()
    gmsh.option.setNumber('Mesh.MeshSizeMin', 300)
    gmsh.option.setNumber('Mesh.MeshSizeMax', 300)
    gmsh.model.mesh.generate(2)

    prior_node_tags, prior_node_coords, _ = gmsh.model.mesh.getNodes()
    prior_coords = prior_node_coords.reshape(-1, 3)
    elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(dim=2)
    tri_idx = list(elem_types).index(2)

    dist_from_pole = np.linalg.norm(prior_coords - np.array([RADIUS, 0, 0]), axis=1)
    metric = np.clip(dist_from_pole * 0.15, 20.0, 300.0)

    gmsh.model.add('planet')
    gmsh.model.setCurrent('planet')
    gmsh.model.occ.addSphere(0, 0, 0, RADIUS)
    gmsh.model.occ.synchronize()

    field_tag = attach_refine_background_field(
        gmsh, prior_node_tags, prior_node_coords,
        elem_tags[tri_idx], elem_node_tags[tri_idx], metric)
    gmsh.model.mesh.field.setAsBackgroundMesh(field_tag)
    gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary', 0)
    gmsh.option.setNumber('Mesh.MeshSizeFromPoints', 0)
    gmsh.option.setNumber('Mesh.MeshSizeFromCurvature', 0)
    gmsh.option.setNumber('Mesh.MeshSizeMin', 1)
    gmsh.option.setNumber('Mesh.MeshSizeMax', 1e6)
    gmsh.model.mesh.generate(2)

    new_node_tags, new_node_coords, _ = gmsh.model.mesh.getNodes()
    print(f'refined mesh: nv={len(new_node_tags)}')
    gmsh.finalize()

'''Generate svFSI ready mesh from a single vtu/vtk file.'''

def separate_surfaces_from_mesh(mesh_file, normals_feature_angle,
                                output_dir=''):
    """
    A VTK script to extract all surfaces which could further be tagged and
    processed in order to isolate surfaces to be used for boundary conditions
    mesh_file: input mesh file (with path and extension)
    """
    import vtk
    import os

    print(mesh_file)
    if mesh_file.endswith('vtu'):
        reader = vtk.vtkXMLUnstructuredGridReader()
    elif mesh_file.endswith('vtk'):
        reader = vtk.vtkUnstructuredGridReader()
    else:
        sys.exit("Only Unstructured Grids Can Be Handled For Now")

    surf_ext = '.vtp'  # output is a vtp file

    reader.SetFileName(mesh_file)
    reader.Update()

    extract_surface = vtk.vtkDataSetSurfaceFilter()
    extract_surface.SetInputConnection(reader.GetOutputPort())
    extract_surface.Update()

    normals = vtk.vtkPolyDataNormals()
    normals.SetFeatureAngle(normals_feature_angle)
    normals.SetInputConnection(extract_surface.GetOutputPort())
    normals.Update()

    connectivity = vtk.vtkConnectivityFilter()
    connectivity.SetInputConnection(normals.GetOutputPort())
    connectivity.ColorRegionsOn()
    connectivity.SetExtractionModeToAllRegions()
    connectivity.Update()
    num_surf = connectivity.GetNumberOfExtractedRegions()

    extracted_regions = connectivity.GetOutput()
    extracted_regions.GetPointData().SetActiveScalars('RegionId')


    if output_dir != '':
        file_dir = output_dir
    else:
        file_dir = os.path.dirname(mesh_file)

    in_file_name = os.path.splitext(os.path.basename(mesh_file))[0]

    output_template = os.path.join(file_dir, in_file_name) + '-surface-{}'

    for i in range(num_surf):

        thresh = vtk.vtkThreshold()
        thresh.SetInputData(extracted_regions)
        thresh.SetInputArrayToProcess(0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS,
                                      vtk.vtkDataSetAttributes.SCALARS, "RegionId")
        thresh.ThresholdBetween(float(i), float(i))
        message = "Lower threshold %d, upper threshold %d" \
                % (thresh.GetLowerThreshold(), thresh.GetUpperThreshold())
        print(message)
        thresh.Update()

        output = thresh.GetOutput()
        output.GetPointData().RemoveArray("Normals")
        output.GetPointData().RemoveArray("RegionId")
        output.GetCellData().RemoveArray("RegionId")
        output.GetCellData().RemoveArray("ModelRegionID")
        FaceID = np.zeros(output.GetNumberOfCells(),dtype=int) + i + 1
        ModelFaceID = vtknp.numpy_to_vtk(FaceID,array_type=vtk.VTK_INT)
        ModelFaceID.SetName("ModelFaceID")
        output.GetCellData().AddArray(ModelFaceID)

        new_surface = vtk.vtkDataSetSurfaceFilter()
        new_surface.SetInputData(output)
        new_surface.Update()

        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetInputData(new_surface.GetOutput())

        writer.SetFileName(output_template.format(i) + surf_ext)
        writer.Update()
        writer.Write()

    return num_surf, output_template, surf_ext

def readvtu(fname):
    # read the ASCII version
    print(fname)
    if fname.endswith('vtu'):
        reader = vtk.vtkXMLUnstructuredGridReader()
    elif fname.endswith('vtk'):
        reader = vtk.vtkUnstructuredGridReader()
    else:
        sys.exit("Only Unstructured Grids Can Be Handled For Now")
    reader.SetFileName(fname)
    reader.Update()
    output = reader.GetOutput()
    return output



import numpy as np 
import vtk
from vtk.util import numpy_support as vtknp
import os

Output1 = readvtu('mesh-complete.mesh.vtu')

GNID = np.linspace(1,Output1.GetNumberOfPoints(),Output1.GetNumberOfPoints(),dtype=int)
GGNID = vtknp.numpy_to_vtk(GNID,array_type=vtk.VTK_INT)
GGNID.SetName("GlobalNodeID")
Output1.GetPointData().AddArray(GGNID)

GEID = np.linspace(1,Output1.GetNumberOfCells(),Output1.GetNumberOfCells(),dtype=int)
GGEID = vtknp.numpy_to_vtk(GEID,array_type=vtk.VTK_INT)
GGEID.SetName("GlobalElementID")
Output1.GetCellData().AddArray(GGEID)

RegionID = np.zeros(Output1.GetNumberOfCells(),dtype=int) + 1
ModelRegionID = vtknp.numpy_to_vtk(RegionID,array_type=vtk.VTK_INT)
ModelRegionID.SetName("ModelRegionID")
Output1.GetCellData().AddArray(ModelRegionID)


# Output
Dir = "ll"
if not os.path.isdir(Dir):
    os.mkdir(Dir)
    os.mkdir(os.path.join(Dir,'mesh-surfaces'))
fname = os.path.join(Dir,'mesh-complete.vtu')

mshwrite = vtk.vtkXMLUnstructuredGridWriter()
mshwrite.SetInputData(Output1)
mshwrite.SetFileName(fname)
mshwrite.Write()

separate_surfaces_from_mesh(fname,50,os.path.join(Dir,'mesh-surfaces'))

# Currently cell position writer is broken, so this is just plotting the vertex positions
# Example written for Nic L as he doesn't have paraview - can extend as desired

import os
import numpy
import vtk
import matplotlib.pyplot as plt
 
output_path = os.path.join(os.getcwd(), 'output')
filename = 'SpheroidTutorial/results_from_time_0/pde_results_oxygen_100.vtu'
 
reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName(os.path.join(output_path, filename))
reader.Update()  # Needed because of GetScalarRange
output = reader.GetOutput()
oxygen = output.GetPointData().GetArray("oxygen")

for i in range(output.GetNumberOfCells()):
   pts = output.GetCell(i).GetPoints()    
   np_pts = numpy.array([pts.GetPoint(i) for i in range(pts.GetNumberOfPoints())])

plt.scatter(np_pts[:,0], np_pts[:,1])
plt.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cell_positions.png'))

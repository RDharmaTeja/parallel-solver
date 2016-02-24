# Visit 2.8.1 log file
ScriptVersion = "2.8.1"
if ScriptVersion != Version():
    print "This script is for VisIt %s. It may not work with version %s" % (ScriptVersion, Version())
visit.ShowAllWindows()
visit.OpenDatabase("./plane-IB.vtk", 0)
# The UpdateDBPluginInfo RPC is not supported in the VisIt module so it will not be logged.
visit.AddPlot("Mesh", "mesh", 1, 1)
visit.DrawPlots()

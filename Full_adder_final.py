#more updates coming soon...
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def run_mumax3(script, name, verbose=False):
    from subprocess import run, PIPE, STDOUT
    from os import path

    scriptfile = name + ".txt"
    outputdir = name + ".out"

    with open(scriptfile, "w") as f:
        f.write(script)

    # call mumax3 to execute this script
    p = run(["mumax3", "-f", scriptfile], stdout=PIPE, stderr=STDOUT)
    if verbose or p.returncode != 0:
        print(p.stdout.decode('UTF-8'))

    if path.exists(outputdir + "/table.txt"):
        table = read_mumax3_table_text(outputdir + "/table.txt")

    else:
        table = None

    #fields = read_mumax3_ovffiles(outputdir)

    return table  #, fields


def read_mumax3_table_text(filename):
    table_file = open(filename, "r")
    tabletext = table_file.read()
    return tabletext


def read_mumax3_table(filename):
    """Puts the mumax3 output table in a pandas dataframe"""

    from pandas import read_table

    table = read_table(filename)
    table.columns = ' '.join(table.columns).split()[1::2]

    return table


def compute_adder(name, frequency, wavelength, eta, thresholdOI, inX, inY, inCi):
    script = f"""
    B0 := 4e-3
    Msat = 1.1e6
    Aex = 18.5e-12
    alpha = 0.004
    Ku1 = 8.3e5
    AnisU = vector(0, 0, 1)
    eta := {eta}
    etaamp := 1/(sqrt(1+eta*eta))
    alpha_0 := 0.004
    alpha_max := 0.75

    theta := 25.84*pi/180
    wavelength := {wavelength}
    frequency := {frequency}
    width := 0.75 * wavelength
    d1 := 6 * wavelength
    d2 := 16 * wavelength
    d3 := 9 * wavelength
    d4 := 1 * wavelength
    d5 := 1 * wavelength
    x1 := 0.5 *d1 *cos(theta)
    h1 := 0.5 *d1*sin(theta)
    x2 := 0.5 *d2 *cos(theta)
    h2 := 0.5 *d2*sin(theta)
    sensorwidth := 0.1 * wavelength
    absorberwidth := 0.75 * wavelength
    deltax := -0.5 * width * sin(theta)
    deltay := -0.5 * width * cos(theta) + 0.5 * width
    
    cellsizex := 2.5e-9
    cellsizey := 2.5e-9
    xmaxsize := 2*x1+d5+d4+10*cellsizex
    ymaxsize := 2*h2 + 10*cellsizey
    setcellsize(cellsizex, cellsizey, 1e-9)
    setgridsize(1536,768,1)

    thresholdOI := {thresholdOI}
    phaseX := 0.0
    phaseY :=0.0
    phaseCi := 0.0
    
    if {inX} == 1{{
        phaseX = pi
    }}
    if {inY} == 1{{
        phaseY = pi
    }}
    if {inCi} == 1{{
        phaseCi = pi
    }}
    phaseI := 0.0

    Region1 := rect(d1,width).RotZ(-theta).transl(-x1-deltax,-h1-deltay,0)
    Region2 := rect(d1,width).RotZ(-theta).transl(x1+deltax,h1+deltay,0)
    Region3 := rect(d1,width).RotZ(theta).transl(x1+deltax,-h1-deltay,0)
    Region4 := rect(d1,width).RotZ(theta).transl(-x1-deltax,h1+deltay,0)
    Region5 := rect(d5,width).transl(-0.5*d5-2*x1-deltax,0,0)
    Region6 := rect(d3,width).transl(0.5*d3 +2*x1+deltax,0,0)
    Region7 := rect(d2,width).RotZ(theta).transl(-x2-deltax,h2+deltay,0)
    Region8 := rect(d2,width).RotZ(-theta).transl(-x2-deltax,-h2-deltay,0)
    Region9 := rect(d2,width).RotZ(-theta).transl(x2+deltax,h2+deltay,0)
    Region10 := rect(d5,width).transl(-0.5*d5 -2*x2-deltax,0,0)
    Region11 := rect(d4,width).transl(0.5*d4 +d3 +2*x1+deltax,0,0)
    InnerWaveguideGeometry := Region1.Add(Region2).Add(Region3).Add(Region4).Add(Region5).Add(Region6)
    OuterWaveguideGeometry := Region7.Add(Region8).Add(Region9).Add(Region10).Add(Region11)

    RegionX := rect(width, width).RotZ(pi/4).transl(0,d1 * sin(theta) - width * sin(theta) + deltay, 0)
    RegionY := rect(width, width).RotZ(pi/4).transl(0,-d1 * sin(theta) + width * sin(theta) - deltay, 0)
    RegionI := rect(width, width).RotZ(pi/4).transl(0,-d2 * sin(theta) + width * sin(theta) - deltay, 0)
    RegionCi:= rect(width, width).RotZ(pi/4).transl(0,d2 * sin(theta) - width * sin(theta) + deltay, 0)
    RegionO := rect(sensorwidth, width).transl(-d5-2*x1-deltax+0.5*sensorwidth, 0,0)
    RegionSum := rect(sensorwidth, width).transl(-d5-2*x2-deltax+0.5*sensorwidth, 0,0)
    RegionCo:= rect(sensorwidth, width).transl(d4+d3+2*x1+deltax-0.5*sensorwidth,0,0)
    Input_geometry := RegionX.Add(RegionY).Add(RegionI).Add(RegionCi)
    Output_geometry:= RegionO.Add(RegionSum).Add(RegionCo)
    
    Absorber1o:= rect(absorberwidth,width).transl(-d5-2*x1-deltax-0.5*absorberwidth,0,0)
    Absorber2o:= rect(absorberwidth,width).transl(-d5-2*x1-deltax-1.5*absorberwidth,0,0)
    Absorber3o:= rect(absorberwidth,width).transl(-d5-2*x1-deltax-2.5*absorberwidth,0,0)
    Absorber4o:= rect(absorberwidth,width).transl(-d5-2*x1-deltax-3.5*absorberwidth,0,0)
    AbsorberO_geom := Absorber1o.Add(Absorber2o).Add(Absorber3o).Add(Absorber4o)

    Absorber1s:= rect(absorberwidth,width).transl(-d5-2*x2-deltax-0.5*absorberwidth,0,0)
    Absorber2s:= rect(absorberwidth,width).transl(-d5-2*x2-deltax-1.5*absorberwidth,0,0)
    Absorber3s:= rect(absorberwidth,width).transl(-d5-2*x2-deltax-2.5*absorberwidth,0,0)
    Absorber4s:= rect(absorberwidth,width).transl(-d5-2*x2-deltax-3.5*absorberwidth,0,0)
    AbsorberS_geom := Absorber1s.Add(Absorber2s).Add(Absorber3s).Add(Absorber4s)
    
    Absorber1co:= rect(absorberwidth,width).transl(d4+d3+2*x1+deltax+0.5*absorberwidth,0,0)
    Absorber2co:= rect(absorberwidth,width).transl(d4+d3+2*x1+deltax+1.5*absorberwidth,0,0)
    Absorber3co:= rect(absorberwidth,width).transl(d4+d3+2*x1+deltax+2.5*absorberwidth,0,0)
    Absorber4co:= rect(absorberwidth,width).transl(d4+d3+2*x1+deltax+3.5*absorberwidth,0,0)
    AbsorberCo_geom := Absorber1co.Add(Absorber2co).Add(Absorber3co).Add(Absorber4co)
    
    Absorber_geom := AbsorberS_geom.Add(AbsorberO_geom).Add(AbsorberCo_geom)
    
    geometry_full_adder := InnerWaveguideGeometry.Add(OuterWaveguideGeometry).Add(Absorber_geom)

    SetGeom(geometry_full_adder)
    saveas(geom, "adder_geom")
    
    DefRegion(8, Absorber1o)
    DefRegion(9, Absorber2o)
    DefRegion(10, Absorber3o)
    DefRegion(11, Absorber4o)
    
    DefRegion(12, Absorber1s)
    DefRegion(13, Absorber2s)
    DefRegion(14, Absorber3s)
    DefRegion(15, Absorber4s)
    
    DefRegion(16, Absorber1co)
    DefRegion(17, Absorber2co)
    DefRegion(18, Absorber3co)
    DefRegion(19, Absorber4co)

    delta_alpha := alpha_max - alpha_0
    alpha.SetRegion(8, alpha_0 + (alpha_max - alpha_0) * Pow(0.25,2))
    alpha.SetRegion(9, alpha_0 + (alpha_max - alpha_0) * Pow(0.5,2))
    alpha.SetRegion(10, alpha_0 + (alpha_max - alpha_0) * Pow(0.75,2))
    alpha.SetRegion(11, alpha_max)
    
    alpha.SetRegion(12, alpha_0 + (alpha_max - alpha_0) * Pow(0.25,2))
    alpha.SetRegion(13, alpha_0 + (alpha_max - alpha_0) * Pow(0.5,2))
    alpha.SetRegion(14, alpha_0 + (alpha_max - alpha_0) * Pow(0.75,2))
    alpha.SetRegion(15, alpha_max)
    
    alpha.SetRegion(16, alpha_0 + (alpha_max - alpha_0) * Pow(0.25,2))
    alpha.SetRegion(17, alpha_0 + (alpha_max - alpha_0) * Pow(0.5,2))
    alpha.SetRegion(18, alpha_0 + (alpha_max - alpha_0) * Pow(0.75,2))
    alpha.SetRegion(19, alpha_max)
    
    m = uniform(0,0,1)
    relax()
    
    saveas(m, "m_start")

    DefRegion(1, RegionX)
    DefRegion(2, RegionY)
    DefRegion(3, RegionCi)
    
    DefRegion(4, RegionO)
    DefRegion(5, RegionI)
    
    DefRegion(6, RegionSum)
    DefRegion(7, RegionCo)
    
    B_ext.SetRegion(1, vector(etaamp * B0 * sin(2*pi*frequency*t+phaseX), -eta * etaamp * B0 * cos(2*pi*frequency*t+phaseX), 0))
    B_ext.SetRegion(2, vector(etaamp * B0 * sin(2*pi*frequency*t+phaseY), -eta * etaamp * B0 * cos(2*pi*frequency*t+phaseY), 0))
    B_ext.SetRegion(3, vector(etaamp * B0 * sin(2*pi*frequency*t+phaseCi), -eta * etaamp * B0 * cos(2*pi*frequency*t+phaseCi), 0))
    
    sensorsize := (sensorwidth * width) / (cellsizex * cellsizey)
    inputsize := (width*width) / (cellsizex * cellsizey)

    mxO := m.Comp(0).Region(4)
    myO := m.Comp(1).Region(4)
    mxS := m.Comp(0).Region(6)
    myS := m.Comp(1).Region(6)
    mxCo:= m.Comp(0).Region(7)
    myCo:= m.Comp(1).Region(7)
    mxX:= m.Comp(0).Region(1)
    myX:= m.Comp(1).Region(1)
    mxY:= m.Comp(0).Region(2)
    myY:= m.Comp(1).Region(2)
    mxCi:= m.Comp(0).Region(3)
    myCi:= m.Comp(1).Region(3)
    
    B_ext.SetRegion(5, vector(B0 * sin(2*pi*frequency*t) * (2*heaviside(sqrt((Sum(Mul(mxO,mxO)) + Sum(Mul(myO,myO)))/sensorsize) - thresholdOI) - 1), -B0 * cos(2*pi*frequency*t) * (2*heaviside(sqrt((Sum(Mul(mxO,mxO)) + Sum(Mul(myO,myO)))/sensorsize) - thresholdOI) - 1),0))

    TableAdd(m.Region(4))
    TableAdd(m.Region(7))
    TableAdd(m.Region(6))
    TableAdd(m.Region(1))
    TableAddVar(sqrt((Sum(Mul(mxO,mxO)) + Sum(Mul(myO,myO)))/sensorsize),"amplitudeO","")
    TableAddVar(sqrt((Sum(Mul(mxS,mxS)) + Sum(Mul(myS,myS)))/sensorsize),"amplitudeS","")
    TableAddVar(sqrt((Sum(Mul(mxX,mxX)) + Sum(Mul(myX,myX)))/inputsize),"amplitudeX","")
    TableAddVar(sqrt((Sum(Mul(mxY,mxY)) + Sum(Mul(myY,myY)))/inputsize),"amplitudeY","")
    TableAddVar(sqrt((Sum(Mul(mxCi,mxCi)) + Sum(Mul(myCi,myCi)))/inputsize),"amplitudeCi","")
TableAddVar(Sum(mxCo)/sensorsize, "mxCo", "")
TableAddVar(Sum(myCo)/sensorsize, "myCo", "")
Tableautosave(1e-12)
//Autosave(m, 1e-12)

run(3e-9)
    """
    table = run_mumax3(script, name, verbose=False)

    out_file = open("./" + name + ".out/data.txt", "w")
    out_file.write(script)

    return table


#(X, Y, Ci)

def write_to_folder(name, folder_name, folder_name_out):
    input_path = "./" + folder_name + ".out" + r"/table.txt"
    input_file = open(input_path, "r")
    table = input_file.read()
    output_path = "./" + folder_name + ".tables/" + folder_name_out + "/" + name + ".txt"
    output_file = open(output_path, "w")
    output_file.write(table)


def write_table_to_folder(table, folder_name, sim_number, table_number, frequency, eta):
    output_path = "./" + folder_name + ".tables/sim" + sim_number + f"/table_{frequency / 1e9:.2f}GHz_eta{eta:.2f}_" + table_number + ".txt"
    output_file = open(output_path, "w")
    output_file.write(table)


def mumax3_output_path(out_folder):
    return "./" + out_folder + ".out/table.txt"


def output_path_text(out_folder, table_number, sim_number):
    path = "./" + out_folder + ".tables/sim" + sim_number + f"/table_{frequency / 1e9:.2f}GHz_eta{eta:.2f}_" + table_number + ".txt"
    return path


def output_path_fig(out_folder, table_number, sim_number, frequency, eta):
    path = "./" + out_folder + ".tables/sim" + sim_number + f"/figure_{frequency / 1e9:.2f}GHz_eta{eta:.2f}_" + table_number + ".png"
    return path


def plot_and_save(out_folder, table_number, sim_number, frequency, eta):
    table_adder = read_mumax3_table(output_path_text(out_folder, table_number, sim_number))

    plt.figure()

    nanosecond = 1e-9
    plt.plot(table_adder["t"] / nanosecond, table_adder["amplitudeO"], label="O")
    plt.plot(table_adder["t"] / nanosecond, table_adder["amplitudeS"], label="Sum")
    plt.plot(table_adder["t"] / nanosecond, table_adder["amplitudeX"], label="X")
    plt.plot(table_adder["t"] / nanosecond, table_adder["amplitudeY"], label="Y")
    plt.plot(table_adder["t"] / nanosecond, table_adder["amplitudeCi"], label="Carry-In")

    plt.xlabel("Time (ns)")
    plt.ylabel("Amplitude")
    plt.legend()

    plt.savefig(output_path_fig(out_folder, table_number, sim_number, frequency, eta), dpi=300, bbox_inches="tight")


def save_data(folder_name, sim_number):
    file = open("./" + folder_name + ".out/data.txt", "r")
    data = file.read()
    out_file = open("./" + folder_name + ".tables/sim" + sim_number + "/data.txt")
    out_file.write(data)


def compute_table(folder_name, sim_number, frequency, wavelength, eta, thresholdOI, inX, inY, inCi):
    table = compute_adder(folder_name, frequency, wavelength, eta, thresholdOI, inX, inY, inCi)
    table_number = str(inX) + str(inY) + str(inCi)
    #table_file = open(mumax3_output_path(folder_name), "r")
    #table = table_file.read()
    write_table_to_folder(table, folder_name, sim_number, table_number, frequency, eta)
    plot_and_save(folder_name, table_number, sim_number, frequency, eta)
    save_data(folder_name, sim_number)


def run_simulation(folder_name, frequency, wavelength, eta, thresholdOI, sim_number):
    compute_table(folder_name, sim_number, frequency, wavelength, eta, thresholdOI, 0, 0, 0)
    compute_table(folder_name, sim_number, frequency, wavelength, eta, thresholdOI, 1, 0, 0)
    compute_table(folder_name, sim_number, frequency, wavelength, eta, thresholdOI, 1, 1, 0)
    compute_table(folder_name, sim_number, frequency, wavelength, eta, thresholdOI, 0, 0, 1)
    compute_table(folder_name, sim_number, frequency, wavelength, eta, thresholdOI, 1, 0, 1)
    compute_table(folder_name, sim_number, frequency, wavelength, eta, thresholdOI, 1, 1, 1)
    sim_number_increment(folder_name)


def sim_number_increment(folder_name):
    sim_file = open("./" + folder_name + ".tables/sim_number_save.txt", "r")
    sim_number_str2 = sim_file.read()
    sim_number_int = int(sim_number_str2)
    out_sim_file = open("./" + folder_name + ".tables/sim_number_save.txt", "w")
    out_sim_file.write(str(sim_number_int + 1))


frequency = 9e9
wavelength = 80.4781e-9
eta = 0.5
thresholdOI = 0.01
sim_file_in = open("./Full_adder_final.tables/sim_number_save.txt", "r")
sim_number_str = sim_file_in.read()
sim_number = int(sim_number_str)

run_simulation("Full_adder_final", frequency, wavelength, eta, thresholdOI, sim_number)

#table000 = compute_adder(folder_name, frequency, wavelength, eta, thresholdOI,0,0,0)
#write_table_to_folder(table000, folder_name, "table000")
#table100 = compute_adder(folder_name, frequency, wavelength, eta, thresholdOI, 1, 0, 0)
#write_table_to_folder(table100, folder_name, "table100")
#table110 = compute_adder(folder_name, frequency, wavelength, eta, thresholdOI, 1, 1, 0)
#write_table_to_folder(table110, folder_name, "table110")
#table001 = compute_adder(folder_name, frequency, wavelength, eta, thresholdOI, 0, 0, 1)
#write_table_to_folder(table001, folder_name, "table001")
#table101 = compute_adder(folder_name, frequency, wavelength, eta, thresholdOI, 1, 0, 1)
#write_table_to_folder(table101, folder_name, "table101")
#table111 = compute_adder(folder_name, frequency, wavelength, eta, thresholdOI, 1, 1, 1)
#write_table_to_folder(table111, folder_name, "table111")


#print(table_adder)

#plt.figure()

#nanosecond = 1e-9
#plt.plot(table_adder["t"] / nanosecond, table_adder["amplitudeO"], label="O")
#plt.plot(table_adder["t"] / nanosecond, table_adder["amplitudeS"], label="Sum")
#plt.plot(table_adder["t"] / nanosecond, table_adder["amplitudeX"], label="X")
#plt.plot(table_adder["t"] / nanosecond, table_adder["amplitudeY"], label="Y")
#plt.plot(table_adder["t"] / nanosecond, table_adder["amplitudeCi"], label="Carry-In")

#plt.xlabel("Time (ns)")
#plt.ylabel("Amplitude")
#plt.legend()

#plt.show()

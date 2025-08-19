<<<<<<< HEAD
import tkinter
import customtkinter
import json
import matplotlib.pyplot as plt
from matplotlib.backend_bases import key_press_handler
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.figure import Figure
from matplotlib.patches import Polygon, FancyArrowPatch
from matplotlib.gridspec import GridSpec
import numpy as np
import datetime
from fipy import Variable, FaceVariable, CellVariable, Grid1D, ImplicitSourceTerm, TransientTerm, DiffusionTerm, Viewer
from fipy.tools import numerix
from scipy.interpolate import CubicSpline

plt.switch_backend('agg')
DIVISION_HORA = 2

# VALIDACIONES DE PARÁMETROS DE ENTRADA    
def validacion_fecha(var):
    try:
        dia = var.timetuple().tm_yday
        mes = var.timetuple().tm_mon
        return True
    except:
        return False
def validacion_zona_horaria(var):
    try:
        if -13 < int(var) < 13:
            return True
        return False
    except:
        return False
def validacion_latitud(var):
    try:
        if -90 <= float(var) <= 90:
            return True
        return False
    except:
        return False
def validacion_orientacion(var):
    try:
        if 0 <= float(var) <= 360:
            return True
        return False
    except:
        return False
def validacion_inclinacion(var):
    try:
        if 0 <= float(var) <= 90:
            return True
        return False
    except:
        return False
def validacion_propiedades_termicas(var):
    try:
        if float(var) > 0:
            return True
        return False
    except:
        return False
def validacion_temperaturas(var):
    try:
        if float(var):
            return True
        return False
    except:
        return False
def validacion_capas(var):
    if var.isdigit():
        return True
    return False
def validacion_propiedades(var):
    if var.isdigit() or var == ".":
        return True
    return False

def solver(espesores, rho_cp, conductividades, difusividades, h_i, h_e, T_e_max, rango_temp, T_i, zona_horaria, fecha, latitud, orientacion, inclinacion):
    """
    Esta función se dedica a la solución numérica del problema
    """

    # MURO
    """
    Se define el tamaño total y el número de puntos a lo largo de este (función de las capas del muro)
    """
    L = sum(espesores) # LARGO DEL MURO
    N_PUNTOS = 1000*int(L/min(espesores)) # NÚMERO DE PUNTOS
    DX = L / N_PUNTOS # DISTANCIA ENTRE PUNTOS
    mesh = Grid1D(nx=N_PUNTOS, dx=DX) # MALLADO DEL MURO

    # TIEMPO
    """
    Se define el intervalo de tiempo y los pasos (1 hora y 24 pasos para considerar un día completo)
    """
    N_DIAS = 5
    DT = 3600/DIVISION_HORA
    PASOS = int(24*3600/DT)*N_DIAS

    # VARIABLES PARA EL MODELO
    """
    Las variables que tendremos en cuenta son el tiempo, la temperatura y la difusividad termica
    """
    tiempo = Variable() # VARIABLE TIEMPO
    temperatura = CellVariable(name="Temperatura", mesh=mesh) # VARIABLE TEMPERATURA
    D = FaceVariable(mesh=mesh) # DIFUSIVIDADES PARA CADA CARA
    K_cell = CellVariable(mesh=mesh) # CONDUCTIVIDADES PARA CADA NODO
    X = mesh.faceCenters[0] # PUNTOS DE LAS CARAS
    RHO_C = CellVariable(mesh=mesh) # PRODUCTO DENSIDAD Y CALOR ESPECIFICO PARA CADA NODO

    for j in range(N_PUNTOS ): 
        for i in range(len(espesores)):
            if X[j] >= sum(espesores[0:i]):
                D[j] = difusividades[i]
                K_cell[j] = conductividades[i]
                RHO_C[j] = conductividades[i]/difusividades[i]

    K = FaceVariable(mesh=mesh) # CONDUCTIVIDADES PARA CADA CARA
    K.setValue(K_cell.harmonicFaceValue)

    # CONDICIONES DE BORDE
    """
    Se tendrán en cuenta efectos convectivos y de radiación, temperatura externa variable y temperatura interior fija
    """
    # TEMPERATURA EXTERNA VARIABLE
    RANGO_TEMPERATURA = rango_temp
    TEMPERTURA_EXTERIOR_MAXIMA = T_e_max
    t_aprox = numerix.linspace(0,24*3600,25)
    PORCENTAJES = numerix.array([0.82, 0.87,0.92,0.96,0.99,1,0.98,0.93,0.84,0.71,0.56,0.39,0.23,0.11,0.03,0,0.03,0.10,0.21,0.34,0.47,0.58,0.68,0.76,0.82])
    cs_aprox = CubicSpline(t_aprox, PORCENTAJES)
    porcentaje = Variable()
    temperatura_exterior = TEMPERTURA_EXTERIOR_MAXIMA - RANGO_TEMPERATURA*porcentaje

    # RADIACION SOLAR
    ORIENTACION = orientacion
    LATITUD = latitud
    INCLINACION = inclinacion
    n_dia = fecha.timetuple().tm_yday
    n_mes = fecha.timetuple().tm_mon
    N_dia = numerix.deg2rad((n_dia-1)*(360/365)) 
    angulo_h = numerix.deg2rad((tiempo % (24*3600) -12*3600)/(24*3600)*360) # ÁNGULO HORARIO (UTC)
    angulo_delta = numerix.deg2rad(0.3963723-22.9132745*numerix.cos(N_dia) + 4.0254304*numerix.sin(N_dia) - 0.3872050*numerix.cos(2*N_dia) + 0.05196728*numerix.sin(2*N_dia) - 0.1545267*numerix.cos(3*N_dia) + 0.08479777*numerix.sin(3*N_dia)) # ÁNGULO DE DECLINACIÓN
    CTE_SOLAR_A = [1202, 1187, 1164, 1130, 1106, 1092, 1093, 1107, 1136, 1166, 1190, 1204]
    CTE_SOLAR_B = [0.141, 0.142, 0.149, 0.164, 0.177, 0.185, 0.186, 0.182, 0.165, 0.152, 0.142, 0.141]
    CTE_SOLAR_C = [0.103, 0.104, 0.109, 0.120, 0.130, 0.137, 0.138, 0.134, 0.121, 0.111, 0.106, 0.103]
    G_T = Variable() # RADIACION TOTAL

    # CONVECCIÓN EXTERIOR + RADIACIÓN EXTERIOR
    mask_e = mesh.facesLeft
    MA = numerix.MA
    tmp = MA.repeat(mesh._faceCenters[..., numerix.NewAxis,:], 2, 1)
    cellToFaceDistanceVectors = tmp - numerix.take(mesh._cellCenters, mesh.faceCellIDs, axis=1)
    tmp = numerix.take(mesh._cellCenters, mesh.faceCellIDs, axis=1)
    tmp = tmp[..., 1,:] - tmp[..., 0,:]
    cellDistanceVectors = MA.filled(MA.where(MA.getmaskarray(tmp), cellToFaceDistanceVectors[:, 0], tmp))
    dPf = FaceVariable(mesh=mesh, value=mesh._faceToCellDistanceRatio * cellDistanceVectors)
    n = mesh.faceNormals
    a_e = h_e*n
    b_e = conductividades[0]
    g_e = h_e*temperatura_exterior + G_T
    coef_robin_e = (mask_e * conductividades[0] * n / (-dPf.dot(a_e) + b_e))

    # TEMPERATURA INTERNA FIJA
    temperatura_interior = T_i # TEMPERATURA INTERNA FIJA

    # CONVECCION INTERIOR
    mask_i = mesh.facesRight
    MA = numerix.MA
    tmp = MA.repeat(mesh._faceCenters[..., numerix.NewAxis,:], 2, 1)
    cellToFaceDistanceVectors = tmp - numerix.take(mesh._cellCenters, mesh.faceCellIDs, axis=1)
    tmp = numerix.take(mesh._cellCenters, mesh.faceCellIDs, axis=1)
    tmp = tmp[..., 1,:] - tmp[..., 0,:]
    cellDistanceVectors = MA.filled(MA.where(MA.getmaskarray(tmp), cellToFaceDistanceVectors[:, 0], tmp))
    dPf = FaceVariable(mesh=mesh, value=mesh._faceToCellDistanceRatio * cellDistanceVectors)
    n = mesh.faceNormals
    a_i = h_i*n
    b_i = conductividades[-1]
    g_i = h_i*temperatura_interior
    coef_robin_i = (mask_i * conductividades[-1] * n / (-dPf.dot(a_i) + b_i))

    # SOLUCIÓN
    """
    La solución considera el término difusivo, el transitorio y las condiciones de borde
    """
    t = numerix.linspace(0, 24*DIVISION_HORA, 24*DIVISION_HORA + 1)
    K.setValue(0., where= mask_e | mask_i)
    eq = (TransientTerm(coeff=RHO_C) == DiffusionTerm(coeff=K) + (coef_robin_e * g_e).divergence - ImplicitSourceTerm(coeff=(coef_robin_e * h_e).divergence) 
        + (coef_robin_i * g_i).divergence - ImplicitSourceTerm(coeff=(coef_robin_i * h_i).divergence)) 

    # VARIABLAS A DEVOLVER
    temperaturas_pasos = []
    temperatura_exterior_pasos = []
    temperaturas_superficie_externa = []
    temperaturas_superficie_interna = []
    temperaturas_superficie_intermedia = []
    flujo_conduccion_entrada = []
    flujo_conduccion_salida = []
    flujo_entrada = []
    flujo_salida = []
    flujo = []

    # LOOP DE TIEMPO
    for paso in range(PASOS + 1):
        tiempo.setValue(paso * DT)  
        porcentaje.setValue(cs_aprox(paso*DT % (24*3600)))
        angulo_beta = numerix.arcsin(numerix.clip((numerix.cos(angulo_delta)*numerix.cos(LATITUD)*numerix.cos(angulo_h)+numerix.sin(angulo_delta)*numerix.sin(LATITUD)), -1, 1)) # ÁNGULO DE ALTITUD SOLAR
        if angulo_h > 0:
            angulo_phi = -numerix.arccos(numerix.clip((numerix.sin(angulo_delta)*numerix.cos(LATITUD)-numerix.sin(LATITUD)*numerix.cos(angulo_delta)*numerix.cos(angulo_h))/(numerix.cos(angulo_beta)), -1, 1)) + 2*numerix.pi
        else:
            angulo_phi = numerix.arccos(numerix.clip((numerix.sin(angulo_delta)*numerix.cos(LATITUD)-numerix.sin(LATITUD)*numerix.cos(angulo_delta)*numerix.cos(angulo_h))/(numerix.cos(angulo_beta)), -1, 1)) # ÁNGULO SOLAR AZIMUTUAL
        angulo_gamma = numerix.absolute(angulo_phi - ORIENTACION) # ÁNGULO SOLAR AZIMUTUAL RESPECTO A LA SUPERFICIE
        angulo_theta = numerix.arccos(numerix.clip(numerix.cos(angulo_beta)*numerix.cos(angulo_gamma)*numerix.sin(INCLINACION) + numerix.sin(angulo_beta)*numerix.cos(INCLINACION), -1, 1)) # ÁNGULO DE INCIDENCIA
        if numerix.sin(angulo_beta) <= .01:
            G_ND = 0
            G_T.setValue(0)
        else:
            G_ND = CTE_SOLAR_A[n_mes-1]/(numerix.exp(CTE_SOLAR_B[n_mes-1]/numerix.sin(angulo_beta)))
            if numerix.rad2deg(INCLINACION) == 90:
                if numerix.cos(angulo_theta) <= - 0.2:
                    G_T.setValue(G_ND*(max(numerix.cos(angulo_theta), 0) + 0.45*CTE_SOLAR_C[n_mes-1]))
                else:
                    G_T.setValue(G_ND*(max(numerix.cos(angulo_theta), 0) + (0.55+0.437*numerix.cos(angulo_theta)+0.313*numerix.cos(angulo_theta)**2)*CTE_SOLAR_C[n_mes-1]))
            else:
                 G_T.setValue(G_ND*(max(numerix.cos(angulo_theta), 0) + (1+numerix.cos(INCLINACION))/2*CTE_SOLAR_C[n_mes-1]))

        eq.solve(var=temperatura, dt=DT)
        temperaturas_pasos.append(temperatura.faceValue.copy())
        temperatura_exterior_pasos.append(temperatura_exterior.copy())
        temperaturas_superficie_externa.append(temperaturas_pasos[paso][0])
        temperaturas_superficie_interna.append(temperaturas_pasos[paso][-1])
        temperaturas_superficie_intermedia.append(temperaturas_pasos[paso][int(len(X)/2)])
        flujo_conduccion_entrada.append(-conductividades[0]*temperatura.faceGrad[0][1].copy())
        flujo_conduccion_salida.append(-conductividades[-1]*temperatura.faceGrad[0][-2].copy())
        flujo_entrada.append(h_e*(temperatura_exterior.copy() - temperaturas_superficie_externa[paso]) + G_T.copy())
        flujo_salida.append(h_i*(-temperatura_interior + temperaturas_superficie_interna[paso]))
        flujo_paso = FaceVariable(mesh=mesh,value=-K*temperatura.faceGrad[0].copy())
        flujo_paso[0] = flujo_paso[1]
        flujo_paso[-1] = flujo_paso[-2]
        flujo.append(flujo_paso) # CORRECCIÓN DEL FLUJO

    # SE RECORTAN LOS VALORES
    temperatura_exterior_pasos = temperatura_exterior_pasos[-24*DIVISION_HORA-1:]
    temperaturas_pasos = temperaturas_pasos[-24*DIVISION_HORA-1:]
    temperaturas_superficie_externa = temperaturas_superficie_externa[-24*DIVISION_HORA-1:]
    temperaturas_superficie_interna = temperaturas_superficie_interna[-24*DIVISION_HORA-1:]
    temperaturas_superficie_intermedia = temperaturas_superficie_intermedia[-24*DIVISION_HORA-1:]
    flujo_conduccion_entrada = flujo_conduccion_entrada[-24*DIVISION_HORA-1:]
    flujo_conduccion_salida = flujo_conduccion_salida[-24*DIVISION_HORA-1:]
    flujo_entrada = flujo_entrada[-24*DIVISION_HORA-1:]
    flujo_salida = flujo_salida[-24*DIVISION_HORA-1:]
    flujo = flujo[-24*DIVISION_HORA-1:]

    
    return [temperaturas_pasos, temperaturas_superficie_externa, temperaturas_superficie_intermedia, temperaturas_superficie_interna, flujo_conduccion_entrada, 
            flujo_conduccion_salida, flujo, X*1000, t/DIVISION_HORA]


# ALERTA
class Alerta(customtkinter.CTkLabel):
    def __init__(self, master, text, type_of_alert,  time):
        super().__init__(master)
        if type_of_alert == 0:
            color = "#007a60"
        else:
            color = "#6e0505"
        self.configure(text=text, text_color="white", fg_color=color)
        self.place(relwidth=1, relx=0, rely=0, anchor="nw", relheight=0.05)
        self.after(time, lambda: self.destroy())

# DATOS DE ENTRADA
class Parametros_entrada:
    def __init__(self):
        
        # CONDICIONES GEOGRÁFICAS
        self.FECHA = 0 # DIA Y MES
        self.ZONA_HORARIA = 0 
        self.LATITUD = 0 
        self.ORIENTACION = 0 # ORIENTACIÓN DEL MURO RESPECTO AL NORTE
        self.INCLINACION = 0 # SUPERFICIE VERTICAL, INCLINADA U HORIZONTAL
        # CONDICIONES INTERNAS
        self.TEMPERATURA_INTERIOR = 0

        # CONDICIONES EXTERNAS
        self.TEMPERATURA_EXTERNA_MAXIMA = 0
        self.RANGO_TEMPERATURA = 0

        # DEFINICIÓN MURO
        self.N_CAPAS = 0
        self.ESPESORES = []
        self.DENSIDADES = []
        self.CALORES_ESPECIFICOS = []
        self.CONDUCTIVIDADES = []
        self.DIFUSIVIDADES = []
        self.PRODUCTO_RHO_CP = []

        # CONVECCIÓN INTERIOR
        self.COEFICIENTE_CONVECCION_INTERIOR = 0

        # CONVECCIÓN EXTERIOR
        self.COEFICIENTE_CONVECCION_EXTERIOR = 0

# PESTAÑA DATOS DE ENTRADA
class Tab_datos_entrada(customtkinter.CTkFrame):
    def __init__(self, master, controller):
        super().__init__(master)
        self.controller = controller
        self.datos_guardados = False
        self.parametros = Parametros_entrada()
        master.grid_columnconfigure(0, weight = 1)
        self.grid_columnconfigure((0,1,2,3,4,5), weight=1)
        
        # FECHA
        self.label_fecha = customtkinter.CTkLabel(self, text="FECHA: ")
        self.label_fecha.grid(row=0, column=0, padx=20, pady=0, sticky="w")
        self.fecha = customtkinter.CTkEntry(self, placeholder_text="DD/MM")
        self.fecha.grid(row=0, column=1, padx=20, pady=0, sticky="ew")

        # ZONA HORARIA
        self.label_zona_horaria = customtkinter.CTkLabel(self, text="ZONA HORARIA: (RESPECTO A UTC)")
        self.label_zona_horaria.grid(row=1, column=0, padx=20, pady=0, sticky="w")
        self.zona_horaria = customtkinter.CTkEntry(self, placeholder_text="-12 a 12")
        self.zona_horaria.grid(row=1, column=1, padx=20, pady=0, sticky="ew")

        # LATITUD
        self.label_latitud = customtkinter.CTkLabel(self, text="LATITUD:")
        self.label_latitud.grid(row=2, column=0, padx=20, pady=0, sticky="w")
        self.latitud = customtkinter.CTkEntry(self, placeholder_text="- 90° a 90°")
        self.latitud.grid(row=2, column=1, padx=20, pady=0, sticky="ew")

        # ORIENTACION
        self.label_orientacion = customtkinter.CTkLabel(self, text="ORIENTACIÓN MURO: ")
        self.label_orientacion.grid(row=3, column=0, padx=20, pady=0, sticky="w")
        self.orientacion = customtkinter.CTkEntry(self, placeholder_text="0° a 360°")
        self.orientacion.grid(row=3, column=1, padx=20, pady=0, sticky="ew")

        # INCLINACION
        self.label_inclinacion = customtkinter.CTkLabel(self, text="INCLINACIÓN MURO: ")
        self.label_inclinacion.grid(row=4, column=0, padx=20, pady=0, sticky="w")
        self.inclinacion = customtkinter.CTkEntry(self, placeholder_text="0° a 90°")
        self.inclinacion.grid(row=4, column=1, padx=20, pady=0, sticky="ew")

        # MAPA
        fig = Figure(facecolor='#2b2b2b')
        ax = fig.add_subplot()
        ax.axis('off')
        fig.tight_layout()
        fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
        img = plt.imread("imagen.png")
        img_altura = 1488
        img_ancho = 2976
        ax.imshow(img, extent=[0,img_ancho,0,img_altura])
        ax.set_xlim(0, img_ancho)
        ax.set_ylim(0, img_altura)
        canvas = FigureCanvasTkAgg(fig, master=self) 
        canvas.draw()
        canvas.get_tk_widget().grid(row=0, column=2, rowspan=6, columnspan=4, sticky="nsew")
        poligono_zona_horaria = None
        vector = None
        muro = None
        # RECARGAR EL MAPA CON LOS DATOS CORRECTOS
        def recargar_mapa(poligono_zona_horaria, vector, muro):
            try:
                for patch in ax.patches[:]:  # Use a copy of the list to avoid modifying it while iterating
                    patch.remove()
                for line in ax.lines[:]:
                    line.remove()
                fecha = datetime.datetime.strptime(self.fecha.get() + "/25", "%d/%m/%y")
                latitud = float(self.latitud.get())
                zona_horaria = int(self.zona_horaria.get())
                orientacion = float(self.orientacion.get())
                inclinacion = float(self.inclinacion.get())
                if not validacion_fecha(fecha):
                    raise ValueError("Fecha no valida")
                if not validacion_latitud(latitud):
                    raise ValueError("Latitud no valida")
                if not validacion_zona_horaria(zona_horaria):
                    raise ValueError("Zona horaria no valida")
                if not validacion_orientacion(orientacion): 
                    raise ValueError("Orientacion no valida")
                if not validacion_inclinacion(inclinacion):
                      raise ValueError("Inclinacion no valida")
                
                if zona_horaria == -12:
                    inicio_longitud = zona_horaria*15
                    fin_longitud = zona_horaria*15 + 7.5
                elif zona_horaria == 12:
                    inicio_longitud = zona_horaria*15 -7.5
                    fin_longitud = zona_horaria*15
                else:
                    inicio_longitud = zona_horaria*15 - 7.5
                    fin_longitud = zona_horaria*15 + 7.5

                # POLIGONO ZONA HORARIA 
                x1, y1 = (inicio_longitud/180*img_ancho/2+img_ancho/2, -90/90*img_altura/2+img_altura/2)
                x2, y2 = (fin_longitud/180*img_ancho/2+img_ancho/2, 90/90*img_altura/2+img_altura/2)
                if poligono_zona_horaria:
                    poligono_zona_horaria.remove()
                poligono_zona_horaria = Polygon(
                    [(x1, y1), (x1, y2), (x2, y2), (x2, y1)], 
                    facecolor="#05646b", alpha=0.3, zorder=2
                )
                ax.add_patch(poligono_zona_horaria)
                
                # LATITUD
                ax.plot([0,img_ancho], [latitud/90*img_altura/2+img_altura/2,latitud/90*img_altura/2+img_altura/2], color="#8ccffc", linewidth=2)

                # MURO Y SU ORIENTACIÓN
                mid_longitud = (inicio_longitud + fin_longitud) / 2
                x0, y0 = (mid_longitud/180*img_ancho/2+img_ancho/2, latitud/90*img_altura/2+img_altura/2)
                angulo_radianes = np.radians(orientacion)
                largo_vector = 150 
                largo_muro = 50
                vector_x = x0 + largo_vector * np.sin(angulo_radianes)
                vector_y = y0 + largo_vector * np.cos(angulo_radianes)
                angulo_radianes_perpendicular = angulo_radianes + np.pi / 2 
                muro_x = x0 + largo_muro * np.sin(angulo_radianes_perpendicular)
                muro_y = y0 + largo_muro * np.cos(angulo_radianes_perpendicular)
                if vector:
                    vector.remove()
                if muro:
                    muro.remove()
                vector = FancyArrowPatch(
                    (x0, y0), (vector_x, vector_y), mutation_scale=15, color="#52bff1", lw=1, zorder=10)
                ax.add_patch(vector)
                muro = FancyArrowPatch(
                    (2*x0 - muro_x, 2*y0 - muro_y), (muro_x, muro_y), arrowstyle="-", mutation_scale=25, color="#1C5B80", lw=4, zorder=11)
                ax.add_patch(muro)
                fig.canvas.draw_idle()

            except ValueError as e:
                Alerta(self.master.master.master, text="Also salió mal, revisar los parámetros de entrada geográficos", type_of_alert=-1, time=1000)
        
        # GUARDADO DEL MAPA
        self.guardado_mapa = customtkinter.CTkButton(self, text="Mostrar en mapa", command=lambda: recargar_mapa(poligono_zona_horaria, vector, muro))
        self.guardado_mapa.grid(row=5, column=0, padx=20, pady=5, sticky="ew", columnspan=2)

        # TEMPERATURA EXTERIOR
        self.label_temperatura_exterior = customtkinter.CTkLabel(self, text="TEMPERATURA EXTERIOR MÁXIMA: ")
        self.label_temperatura_exterior.grid(row=6, column=0, padx=20, pady=10, sticky="w")
        self.temperatura_exterior_maxima = customtkinter.CTkEntry(self, placeholder_text="°C")
        self.temperatura_exterior_maxima.grid(row=6, column=1, padx=20, pady=10, sticky="ew")
        self.label_rango_temperatura = customtkinter.CTkLabel(self, text="RANGO DE TEMPERATURA EXTERIOR: ")
        self.label_rango_temperatura.grid(row=6, column=2, padx=20, pady=10, sticky="w")
        self.rango_temperatura = customtkinter.CTkEntry(self, placeholder_text="°C")
        self.rango_temperatura.grid(row=6, column=3, padx=20, pady=10, sticky="ew")

        # CONVECCION EXTERIOR
        self.label_conveccion_exterior = customtkinter.CTkLabel(self, text="COEFICIENTE CONVECCION EXTERIOR: ")
        self.label_conveccion_exterior.grid(row=6, column=4, padx=20, pady=10, sticky="w")
        self.coeficiente_conveccion_exterior = customtkinter.CTkEntry(self, placeholder_text="W/m²K")
        self.coeficiente_conveccion_exterior.grid(row=6, column=5, padx=20, pady=10, sticky="ew")

        # TEMPERATURA INTERIOR
        self.label_temperatura_interior = customtkinter.CTkLabel(self, text="TEMPERATURA INTERIOR: ")
        self.label_temperatura_interior.grid(row=7, column=0, padx=20, pady=10, sticky="w")
        self.temperatura_interior = customtkinter.CTkEntry(self, placeholder_text="°C")
        self.temperatura_interior.grid(row=7, column=1, padx=20, pady=10, sticky="ew")

        # CONVECCION INTERIOR
        self.label_conveccion_interior = customtkinter.CTkLabel(self, text="COEFICIENTE CONVECCIÓN INTERIOR: ")
        self.label_conveccion_interior.grid(row=7, column=2, padx=20, pady=10, sticky="w")
        self.coeficiente_conveccion_interior = customtkinter.CTkEntry(self, placeholder_text="W/m²K")
        self.coeficiente_conveccion_interior.grid(row=7, column=3, padx=20, pady=10, sticky="ew")

        # NÚMERO DE CAPAS DEL MURO
        self.label_n_capas = customtkinter.CTkLabel(self, text="NÚMERO DE CAPAS: ")
        self.label_n_capas.grid(row=7, column=4, padx=20, pady=10, sticky="w")
        self.n_capas = customtkinter.CTkEntry(self, placeholder_text="", validate="key", validatecommand=(self.register(validacion_capas), "%S"))
        self.n_capas.grid(row=7, column=5, padx=20, pady=10, sticky="ew")
        self.propiedades_muro = ["N° DE CAPA", "ESPESOR [mm]", "DENSIDAD [kg/m³]", "CALOR ESPECÍFICO [J/kgK]", "CONDUCTIVIDAD TÉRMICA [W/mK]"]
        self.labels_temporales = []
        self.labels_capas = []
        self.espesores_capas = []
        self.densidades_capas = []
        self.calores_especificos_capas = []
        self.conductividades_termicas_capas = []

        self.boton_n_capas = customtkinter.CTkButton(self, text="Seleccionar capas", command=self.cargar_capas)
        self.boton_n_capas.grid(row=8, column=0, padx=20, pady=10, sticky="ew", columnspan=6)

        self.guardar_todo = customtkinter.CTkButton(self, text="Guardar todo", fg_color="#23a15c", hover_color="#136337", command=self.carga_de_datos)

    # TABLA CON DATOS DEL MURO
    def cargar_capas(self):
        if self.n_capas.get() == "":
            return True
        else:
            n_capas = max(1, min(int(self.n_capas.get()), 5))
            if len(self.labels_temporales) == 0:
                for i in range(len(self.propiedades_muro)):
                    tmp = customtkinter.CTkLabel(self, text=self.propiedades_muro[i])
                    tmp.grid(row = 9 , column=i, padx=20, pady=10)
                    self.labels_temporales.append(tmp)
            else:
                for i in range(len(self.labels_capas)):
                    self.labels_capas[i].destroy()
                    self.espesores_capas[i].destroy()
                    self.densidades_capas[i].destroy()
                    self.calores_especificos_capas[i].destroy()
                    self.conductividades_termicas_capas[i].destroy()
                self.labels_capas = []
                self.espesores_capas = []
                self.densidades_capas = []
                self.calores_especificos_capas = []
                self.conductividades_termicas_capas = []
                self.guardar_todo.grid_forget()
            for i in range(n_capas):
                    label_capa = customtkinter.CTkLabel(self, text=i + 1)
                    label_capa.grid(row = 10 + i, column= 0, padx=20, pady=10)
                    self.labels_capas.append(label_capa)
                    label_espesor = customtkinter.CTkEntry(self, validate="key", validatecommand=(self.register(validacion_propiedades), "%S"))
                    label_espesor.grid(row = 10 + i, column = 1, padx=20, pady=10)
                    self.espesores_capas.append(label_espesor)
                    label_densidad = customtkinter.CTkEntry(self, validate="key", validatecommand=(self.register(validacion_propiedades), "%S"))
                    label_densidad.grid(row = 10 + i, column = 2, padx=20, pady=10)
                    self.densidades_capas.append(label_densidad)
                    label_calor_especifico = customtkinter.CTkEntry(self, validate="key", validatecommand=(self.register(validacion_propiedades), "%S"))
                    label_calor_especifico.grid(row = 10 + i, column = 3, padx=20, pady=10)
                    self.calores_especificos_capas.append(label_calor_especifico)
                    label_conductividad = customtkinter.CTkEntry(self, validate="key", validatecommand=(self.register(validacion_propiedades), "%S"))
                    label_conductividad.grid(row = 10 + i, column = 4, padx=20, pady=10)
                    self.conductividades_termicas_capas.append(label_conductividad)
            self.guardar_todo.grid(row=10, column=5, padx=20, pady=10, rowspan=n_capas, sticky="nsew")

    # CARGAR LOS DATOS Y VALIDAR LAS ENTRADAS
    def carga_de_datos(self):
        try:
            self.parametros.FECHA = datetime.datetime.strptime(self.fecha.get() + "/25", "%d/%m/%y")
            self.parametros.LATITUD = numerix.deg2rad(float(self.latitud.get()))
            self.parametros.ZONA_HORARIA = numerix.deg2rad(int(self.zona_horaria.get()))
            self.parametros.ORIENTACION = numerix.deg2rad(float(self.orientacion.get()))
            self.parametros.INCLINACION = numerix.deg2rad(float(self.inclinacion.get()))
            if (not validacion_fecha(self.parametros.FECHA) or not validacion_latitud(self.parametros.LATITUD) or not validacion_zona_horaria(self.parametros.ZONA_HORARIA) or
                not validacion_orientacion(self.parametros.ORIENTACION) or not validacion_inclinacion(self.parametros.INCLINACION)):
                raise ValueError("Datos geográficos inválidos")
            self.parametros.COEFICIENTE_CONVECCION_INTERIOR = float(self.coeficiente_conveccion_interior.get())
            self.parametros.COEFICIENTE_CONVECCION_EXTERIOR = float(self.coeficiente_conveccion_exterior.get())
            self.parametros.TEMPERATURA_INTERIOR = float(self.temperatura_interior.get())
            self.parametros.TEMPERATURA_EXTERNA_MAXIMA = float(self.temperatura_exterior_maxima.get())
            self.parametros.RANGO_TEMPERATURA = float(self.rango_temperatura.get())
            self.parametros.N_CAPAS = max(1, min(int(self.n_capas.get()), 5))
            if (not validacion_propiedades_termicas(self.parametros.COEFICIENTE_CONVECCION_INTERIOR) or not validacion_propiedades_termicas(self.parametros.COEFICIENTE_CONVECCION_EXTERIOR)
            or not validacion_temperaturas(self.parametros.TEMPERATURA_INTERIOR) or not validacion_temperaturas(self.parametros.TEMPERATURA_EXTERNA_MAXIMA) or
            not validacion_temperaturas(self.parametros.RANGO_TEMPERATURA)):
                raise ValueError("Datos térmicos inválidos")
            if len(self.parametros.ESPESORES) != 0:
                self.parametros.ESPESORES = []
                self.parametros.DENSIDADES = []
                self.parametros.CALORES_ESPECIFICOS = []
                self.parametros.CONDUCTIVIDADES = []
                self.parametros.DIFUSIVIDADES = []
                self.parametros.PRODUCTO_RHO_CP = []
            for i in range(self.parametros.N_CAPAS):
                self.parametros.ESPESORES.append(float(self.espesores_capas[i].get())/1000)
                self.parametros.DENSIDADES.append(float(self.densidades_capas[i].get()))
                self.parametros.CALORES_ESPECIFICOS.append(float(self.calores_especificos_capas[i].get()))
                self.parametros.CONDUCTIVIDADES.append(float(self.conductividades_termicas_capas[i].get()))
                self.parametros.DIFUSIVIDADES.append(self.parametros.CONDUCTIVIDADES[i]/(self.parametros.CALORES_ESPECIFICOS[i]*self.parametros.DENSIDADES[i]))
                self.parametros.PRODUCTO_RHO_CP.append(self.parametros.CALORES_ESPECIFICOS[i]*self.parametros.DENSIDADES[i])
                if (not validacion_propiedades_termicas(self.parametros.ESPESORES[i]) or not validacion_propiedades_termicas(self.parametros.DENSIDADES[i])
                or not validacion_propiedades_termicas(self.parametros.CALORES_ESPECIFICOS[i]) or not validacion_propiedades_termicas(self.parametros.CONDUCTIVIDADES[i])):
                    raise ValueError("Datos del muro inválidos")
            self.datos_guardados = True
            self.controller.tab_perfiles_temperatura.boton_simulacion.configure(state=customtkinter.NORMAL)
            self.controller.tab_resultados_termicos.boton_simulacion.configure(state=customtkinter.NORMAL)
            Alerta(self.master.master.master, text="Datos guardados correctamente", type_of_alert=0, time=1000)
        except ValueError as e:
            if str(e) == "Datos geográficos inválidos" or str(e) == "Datos térmicos inválidos":
                Alerta(self.master.master.master, text=str(e), type_of_alert=-1, time=1000)
            else:
                Alerta(self.master.master.master, text="Algo salio mal, revisar los parámetros de entrada", type_of_alert=-1, time=1000)
        except:
            Alerta(self.master.master.master, text="Algo salio mal, revisar los parámetros de entrada", type_of_alert = -1, time=1000)
    
            
# PESTAÑA PERFILES DE TEMPERATURA
class Tab_perfiles_de_temperaturas(customtkinter.CTkFrame):
    def __init__(self, master, controller):
        super().__init__(master)
        self.controller = controller
        master.grid_columnconfigure(0, weight = 1)
        master.grid_rowconfigure(0, weight = 1)
        self.grid_columnconfigure((0,1,2,3,4,5), weight=1)
        self.grid_rowconfigure((0,1,2,3,4,5), weight=1)

        # GRÁFICOS
        plt.style.use("dark_background")
        self.fig = Figure(facecolor='#2b2b2b')
        gs = GridSpec(2, 1, wspace=0.2, hspace=0.4, left=0.1, right=0.9, top=0.95, bottom=0.06)

        # GRÁFICO TEMPERATURA DEL MURO
        self.ax_temp_muro = self.fig.add_subplot(gs[0, 0]) 
        self.ax_temp_muro.set_title("Perfil de temperatura en el muro")
        self.ax_temp_muro.set_xlabel(r"$e$ [mm]")
        self.ax_temp_muro.set_ylabel(r"$T$ [°C]")

        # GRÁFICO TEMPERATURAS EN SUPERFICIES
        self.ax_temp_superficie = self.fig.add_subplot(gs[1, 0])
        self.ax_temp_superficie.set_xticks([0,2,4,6,8,10,12,14,16,18,20,22,24])
        self.ax_temp_superficie.set_title("Temperatura en superficie exterior, interior e intermedia")
        self.ax_temp_superficie.set_xlim(0,24)
        self.ax_temp_superficie.set_xlabel(r"Hora [h]")
        self.ax_temp_superficie.set_ylabel(r"$T$ [°C]")
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=self) 
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=0, column=0, rowspan=5, columnspan=6, sticky="nsew")
        
        # CONTROLES
        self.label_slider = customtkinter.CTkLabel(self, text="HORA A ANALIZAR: 00:00")
        self.label_slider.grid(row=5,column=0, sticky="e", padx=20, pady=10)
        self.slider_simulacion = customtkinter.CTkSlider(self, from_=0, to=24, number_of_steps=24, state=customtkinter.DISABLED, command=self.seleccionar_hora)
        self.slider_simulacion.set(0)
        self.slider_simulacion.grid(row=5,column=1, columnspan=2, sticky="ew", padx=20, pady=10)
        self.boton_simulacion = customtkinter.CTkButton(self, text="Correr simulación", command=self.controller.solucion_numerica, state=customtkinter.DISABLED)
        self.boton_simulacion.grid(row=5, column=4, sticky="ew", padx=20, pady=10)

        self.min_temp = 0
        self.max_temp = 0

    def seleccionar_hora(self, value):
        value = int(value)
        for line in self.ax_temp_muro.lines[:]:
            line.remove()
        self.ax_temp_muro.plot(self.controller.resultados_simulacion[-2], self.controller.resultados_simulacion[0][value*DIVISION_HORA], color="C1")
        if int(value) == 0:
            value = "00"
        self.label_slider.configure(text=f"HORA A ANALIZAR: {int(value)}:00")
        acum_espesor = 0
        for i in self.controller.tab_datos_entrada.parametros.ESPESORES:
            acum_espesor += i
            self.ax_temp_muro.plot([(acum_espesor)*1000,(acum_espesor)*1000], [self.min_temp, self.max_temp], linestyle="dashed", lw=1, color="C3")
        self.canvas.draw_idle()

# PESTAÑA RESULTADOS TÉRMICOS
class Tab_resultados_termicos(customtkinter.CTkFrame):
    def __init__(self, master, controller):
        super().__init__(master)
        self.controller = controller
        master.grid_columnconfigure(0, weight = 1)
        master.grid_rowconfigure(0, weight = 1)
        self.grid_columnconfigure((0,1,2,3,4,5), weight=1)
        self.grid_rowconfigure((0,1,2,3,4,5), weight=1)
        
        # GRÁFICOS
        plt.style.use("dark_background")
        self.fig = Figure(facecolor='#2b2b2b')
        gs = GridSpec(2, 1, wspace=0.2, hspace=0.4, left=0.1, right=0.9, top=0.95, bottom=0.06)

        # GRÁFICO FLUJO DE CALOR EN EL MURO
        self.ax_flujo_calor = self.fig.add_subplot(gs[0, 0]) 
        self.ax_flujo_calor.set_title("Flujo de calor a lo largo del muro")
        self.ax_flujo_calor.set_xlabel(r"$e$ [mm]")
        self.ax_flujo_calor.set_ylabel(r"$Q$ $[W/m^2]$")

        # GRÁFICO GANANCIA y CARGA TÉRMICA
        self.ax_ganancia_carga_termica = self.fig.add_subplot(gs[1, 0])
        self.ax_ganancia_carga_termica.set_title("Ganancia y carga térmica")
        self.ax_ganancia_carga_termica.set_xlim(0,24)
        self.ax_ganancia_carga_termica.set_xticks([0,2,4,6,8,10,12,14,16,18,20,22,24])
        self.ax_ganancia_carga_termica.set_xlabel(r"Hora [h]")
        self.ax_ganancia_carga_termica.set_ylabel(r"$Q$ $[W/m^2]$")
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=self) 
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=0, column=0, rowspan=5, columnspan=6, sticky="nsew")

        # CONTROLES
        self.label_slider = customtkinter.CTkLabel(self, text="HORA A ANALIZAR: 00:00")
        self.label_slider.grid(row=5,column=0, sticky="e", padx=20, pady=10)
        self.slider_simulacion = customtkinter.CTkSlider(self, from_=0, to=24, number_of_steps=24, state=customtkinter.DISABLED, command=self.seleccionar_hora)
        self.slider_simulacion.set(0)
        self.slider_simulacion.grid(row=5,column=1, columnspan=2, sticky="ew", padx=20, pady=10)
        self.boton_simulacion = customtkinter.CTkButton(self, text="Correr simulación", command=self.controller.solucion_numerica, state=customtkinter.DISABLED)
        self.boton_simulacion.grid(row=5, column=4, sticky="ew", padx=20, pady=10)

        self.min_flujo = 0
        self.max_flujo = 0

    def seleccionar_hora(self, value):
        value = int(value)
        for line in self.ax_flujo_calor.lines[:]:
            line.remove()
        self.ax_flujo_calor.plot(self.controller.resultados_simulacion[-2], self.controller.resultados_simulacion[-3][value*DIVISION_HORA], color="C1")
        if int(value) == 0:
            value = "00"
        self.label_slider.configure(text=f"HORA A ANALIZAR: {int(value)}:00")
        acum_espesor = 0
        for i in self.controller.tab_datos_entrada.parametros.ESPESORES:
            acum_espesor += i
            self.ax_flujo_calor.plot([(acum_espesor)*1000,(acum_espesor)*1000],[self.min_flujo, self.max_flujo], linestyle="dashed", lw=1, color="C3")
        self.canvas.draw_idle()

# PESTAÑA COMPARACIÓN MÚLTIPLE
class Tab_comparacion_multiple(customtkinter.CTkFrame):
    def __init__(self, master, controller):
        super().__init__(master)
        self.controller = controller
        master.grid_columnconfigure(0, weight = 1)
        master.grid_rowconfigure(0, weight = 1)
        self.grid_columnconfigure((0,1,2,3,4,5,6), weight=1)
        self.grid_rowconfigure((0,1,2,3,4,5), weight=1)
        self.muros = []

        # BOTÓN CARGA ARCHIVO
        self.boton_carga_archivo = customtkinter.CTkButton(self, text="Cargar archivo .json", fg_color="#23a15c", hover_color="#136337", command=self.cargar_json)
        self.boton_carga_archivo.grid(row=0, column=0, padx=10, sticky="ew")

        # TEXTBOX CON DATOS DE LOS MUROS
        self.texto_json = customtkinter.CTkTextbox(self, width=550, height=300)
        self.texto_json.grid(row=1, column=0, padx=10, sticky="nsew", rowspan=4)

        # CORRER SIMULACION
        self.boton_simulacion = customtkinter.CTkButton(self, text="Correr simulación", command=self.controller.solucion_numerica_multiple, state=customtkinter.DISABLED)
        self.boton_simulacion.grid(row=5, column=0, sticky="ew", padx=10)
        
        # GRÁFICOS
        plt.style.use("dark_background")
        self.fig = Figure(facecolor='#2b2b2b')
        gs = GridSpec(2, 1, wspace=0.2, hspace=0.4, left=0.1, right=0.9, top=0.95, bottom=0.06)


        # GRÁFICO GANANCIA TÉRMICA
        self.ax_ganancia_termica = self.fig.add_subplot(gs[0, 0])
        self.ax_ganancia_termica.set_title("Ganancia térmica")
        self.ax_ganancia_termica.set_xlim(0,24)
        self.ax_ganancia_termica.set_xticks([0,2,4,6,8,10,12,14,16,18,20,22,24])
        self.ax_ganancia_termica.set_xlabel(r"Hora [h]")
        self.ax_ganancia_termica.set_ylabel(r"$Q$ $[W/m^2]$")

        # GRÁFICO CARGA TÉRMICA
        self.ax_carga_termica = self.fig.add_subplot(gs[1, 0])
        self.ax_carga_termica.set_title("Carga térmica")
        self.ax_carga_termica.set_xlim(0,24)
        self.ax_carga_termica.set_xticks([0,2,4,6,8,10,12,14,16,18,20,22,24])
        self.ax_carga_termica.set_xlabel(r"Hora [h]")
        self.ax_carga_termica.set_ylabel(r"$Q$ $[W/m^2]$")
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=self) 
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=0, column=1, rowspan=6, columnspan=5, sticky="nsew")

    # CARGA Y PROCESAMIENTO DEL ARCHIVO .JSON
    def cargar_json(self):
        archivo = tkinter.filedialog.askopenfilename(
            title="Seleccionar archivo JSON",
            filetypes=[("Archivos JSON", "*.json")]
        )
        self.muros = []
        if archivo:
            try:
                with open(archivo, "r", encoding="utf-8") as f:
                    datos = json.load(f)
                    for muro in datos:
                        id_muro = muro["id"]
                        nombre_muro = muro["nombre"]
                        parametros = Parametros_entrada()
                        parametros.FECHA = datetime.datetime.strptime(muro["fecha"] + "/25", "%d/%m/%y")
                        parametros.ZONA_HORARIA = numerix.deg2rad(muro["zona_horaria"])
                        parametros.LATITUD = numerix.deg2rad(muro["latitud"])
                        parametros.ORIENTACION = numerix.deg2rad(muro["orientacion"])
                        parametros.INCLINACION = numerix.deg2rad(muro["inclinacion"])
                        if (not validacion_fecha(parametros.FECHA) or not validacion_latitud(parametros.LATITUD) or not validacion_zona_horaria(parametros.ZONA_HORARIA) or
                            not validacion_orientacion(parametros.ORIENTACION) or not validacion_inclinacion(parametros.INCLINACION)):
                            raise ValueError("Datos geográficos inválidos")
                        parametros.COEFICIENTE_CONVECCION_INTERIOR = muro["coeficiente_conveccion_interior"]
                        parametros.COEFICIENTE_CONVECCION_EXTERIOR = muro["coeficiente_conveccion_exterior"]
                        parametros.TEMPERATURA_INTERIOR = muro["temperatura_interior"]
                        parametros.TEMPERATURA_EXTERNA_MAXIMA = muro["temperatura_exterior"]
                        parametros.RANGO_TEMPERATURA = muro["rango_temperatura_exterior"]
                        parametros.N_CAPAS = muro["n_capas"]
                        
                        if (not validacion_propiedades_termicas(parametros.COEFICIENTE_CONVECCION_INTERIOR) or not validacion_propiedades_termicas(parametros.COEFICIENTE_CONVECCION_EXTERIOR)
                        or not validacion_temperaturas(parametros.TEMPERATURA_INTERIOR) or not validacion_temperaturas(parametros.TEMPERATURA_EXTERNA_MAXIMA) or
                        not validacion_temperaturas(parametros.RANGO_TEMPERATURA)):
                            raise ValueError("Datos térmicos inválidos")
                        for capa in muro["capas"]:
                            parametros.ESPESORES.append(capa["espesor"]/1000)
                            parametros.DENSIDADES.append(capa["densidad"])
                            parametros.CALORES_ESPECIFICOS.append(capa["calor_especifico"])
                            parametros.CONDUCTIVIDADES.append(capa["conductividad"])
                            parametros.DIFUSIVIDADES.append(capa["conductividad"]/(capa["densidad"]*capa["calor_especifico"]))
                            parametros.PRODUCTO_RHO_CP.append(capa["densidad"]*capa["calor_especifico"])
                            if (not validacion_propiedades_termicas(parametros.ESPESORES[-1]) or not validacion_propiedades_termicas(parametros.DENSIDADES[-1])
                            or not validacion_propiedades_termicas(parametros.CALORES_ESPECIFICOS[-1]) or not validacion_propiedades_termicas(parametros.CONDUCTIVIDADES[-1])):
                                raise ValueError("Datos del muro inválidos")
                        self.muros.append([id_muro, nombre_muro, parametros])
                    contenido = json.dumps(datos, indent=4, ensure_ascii=False)
                    self.texto_json.delete("0.0", "end")
                    self.texto_json.insert("0.0", contenido)
                    self.boton_simulacion.configure(state=customtkinter.NORMAL)
                    Alerta(self.master.master.master, text="Archivo .json cargado correctamente", type_of_alert=0, time=1000)
            except Exception as e:
                Alerta(self.master.master.master, text="Error al cargar el archivo .json, revisar los datos ingresados", type_of_alert=-1, time=2000)

# PESTAÑAS
class MyTabView(customtkinter.CTkTabview):
    def __init__(self, master, **kwargs):
        super().__init__(master, **kwargs)
        self.resultados_simulacion = []
        self.resultados_simulacion_multiple = []

        self.grid_columnconfigure(0, weight=1)
        # PESTAÑAS
        self.add("Datos de entrada")
        self.add("Resultados temperatura")
        self.add("Resultados térmicos")
        self.add("Comparación múltiple")

        # WIDGETS DATOS DE ENTRADA
        self.tab_datos_entrada = Tab_datos_entrada(master=self.tab("Datos de entrada"), controller=self)
        self.tab_datos_entrada.grid(row=0, column=0, padx=20, pady=10, sticky="nsew")

        # WIDGETS PERFILES DE TEMPERATURA
        self.tab_perfiles_temperatura = Tab_perfiles_de_temperaturas(master=self.tab("Resultados temperatura"), controller=self)
        self.tab_perfiles_temperatura.grid(row=0, column=0, padx=20, pady=10, sticky="nsew")

        # WIDGETS RESULTADOS TÉRMICOS
        self.tab_resultados_termicos = Tab_resultados_termicos(master=self.tab("Resultados térmicos"), controller=self)
        self.tab_resultados_termicos.grid(row=0, column=0, padx=20, pady=10, sticky="nsew")

        # WIDGETS COMPARACION MULTIPLE
        self.tab_comparacion_multiple = Tab_comparacion_multiple(master=self.tab("Comparación múltiple"), controller=self)
        self.tab_comparacion_multiple.grid(row=0, column=0, padx=20, pady=10, sticky="nsew")

    def solucion_numerica(self):
        try:
            # OBTENER RESULTADOS DE SIMULACIÓN
            self.resultados_simulacion = solver(espesores=self.tab_datos_entrada.parametros.ESPESORES, rho_cp=self.tab_datos_entrada.parametros.PRODUCTO_RHO_CP, 
                        conductividades=self.tab_datos_entrada.parametros.CONDUCTIVIDADES, difusividades=self.tab_datos_entrada.parametros.DIFUSIVIDADES,
                        h_i=self.tab_datos_entrada.parametros.COEFICIENTE_CONVECCION_INTERIOR, h_e=self.tab_datos_entrada.parametros.COEFICIENTE_CONVECCION_EXTERIOR,
                        T_e_max=self.tab_datos_entrada.parametros.TEMPERATURA_EXTERNA_MAXIMA, rango_temp=self.tab_datos_entrada.parametros.RANGO_TEMPERATURA,
                        T_i=self.tab_datos_entrada.parametros.TEMPERATURA_INTERIOR, zona_horaria=self.tab_datos_entrada.parametros.ZONA_HORARIA,
                        fecha=self.tab_datos_entrada.parametros.FECHA, latitud=self.tab_datos_entrada.parametros.LATITUD, 
                        orientacion=self.tab_datos_entrada.parametros.ORIENTACION, inclinacion=self.tab_datos_entrada.parametros.INCLINACION)

            # PERFIL DE TEMPERATURAS
            for line in self.tab_perfiles_temperatura.ax_temp_muro.lines[:]:
                line.remove()
            for line in self.tab_perfiles_temperatura.ax_temp_superficie.lines[:]:
                line.remove()
            
            # PARÁMETROS
            slider_temp_valor = int(self.tab_perfiles_temperatura.slider_simulacion.get())
            self.tab_perfiles_temperatura.min_temp = int(min(min(i) for i in self.resultados_simulacion[0]))
            self.tab_perfiles_temperatura.max_temp = int(max(max(i) for i in self.resultados_simulacion[0])) + 1
           
            # GRÁFICO DE LAS SUPERFICIES
            self.tab_perfiles_temperatura.ax_temp_superficie.plot(self.resultados_simulacion[-1], self.resultados_simulacion[1], color="C0", label="Superficie exterior")
            self.tab_perfiles_temperatura.ax_temp_superficie.plot(self.resultados_simulacion[-1], self.resultados_simulacion[2], color="C2", label="Superficie intermedia")
            self.tab_perfiles_temperatura.ax_temp_superficie.plot(self.resultados_simulacion[-1], self.resultados_simulacion[3], color="C3", label="Superficie interior")
            self.tab_perfiles_temperatura.ax_temp_superficie.legend()

            # GRÁFICO PERFIL DEL MURO
            self.tab_perfiles_temperatura.ax_temp_muro.plot(self.resultados_simulacion[-2], self.resultados_simulacion[0][slider_temp_valor], color="C1")
            self.tab_perfiles_temperatura.ax_temp_muro.set_ylim(self.tab_perfiles_temperatura.min_temp, self.tab_perfiles_temperatura.max_temp)
            self.tab_perfiles_temperatura.ax_temp_muro.set_yticks(np.linspace(self.tab_perfiles_temperatura.min_temp, self.tab_perfiles_temperatura.max_temp, 5))
            self.tab_perfiles_temperatura.ax_temp_muro.set_xlim(self.resultados_simulacion[-2][0],self.resultados_simulacion[-2][-1])
            acum_espesor = 0
            for i in self.tab_datos_entrada.parametros.ESPESORES:
                acum_espesor += i
                self.tab_perfiles_temperatura.ax_temp_muro.plot([(acum_espesor)*1000,(acum_espesor)*1000], [self.tab_perfiles_temperatura.min_temp, self.tab_perfiles_temperatura.max_temp], linestyle="dashed", lw=1, color="C3")
            
            
            self.tab_perfiles_temperatura.slider_simulacion.configure(state=customtkinter.NORMAL)
            self.tab_perfiles_temperatura.canvas.draw_idle()

            # FLUJOS DE CALOR
            for line in self.tab_resultados_termicos.ax_flujo_calor.lines[:]:
                line.remove()
            for line in self.tab_resultados_termicos.ax_ganancia_carga_termica.lines[:]:
                line.remove()

            # PARÁMETROS
            slider_flujo_valor = int(self.tab_resultados_termicos.slider_simulacion.get())
            self.tab_resultados_termicos.min_flujo = int(min(min(i) for i in self.resultados_simulacion[-3]))
            self.tab_resultados_termicos.max_flujo = int(max(max(i) for i in self.resultados_simulacion[-3])) + 1

            # GRÁFICO DE CARGA Y GANANCIA
            self.tab_resultados_termicos.ax_ganancia_carga_termica.plot(self.resultados_simulacion[-1], self.resultados_simulacion[-5], color="C0", label="Ganancia térmica")
            self.tab_resultados_termicos.ax_ganancia_carga_termica.plot(self.resultados_simulacion[-1], self.resultados_simulacion[-4], color="C2", label="Carga térmica")
            self.tab_resultados_termicos.ax_ganancia_carga_termica.legend()
            
            # GRÁFICO FLUJO CALOR
            self.tab_resultados_termicos.ax_flujo_calor.plot(self.resultados_simulacion[-2], self.resultados_simulacion[-3][slider_flujo_valor], color="C1")
            self.tab_resultados_termicos.ax_flujo_calor.set_ylim(self.tab_resultados_termicos.min_flujo, self.tab_resultados_termicos.max_flujo)
            self.tab_resultados_termicos.ax_flujo_calor.set_yticks(np.linspace(self.tab_resultados_termicos.min_flujo, self.tab_resultados_termicos.max_flujo, 5))
            self.tab_resultados_termicos.ax_flujo_calor.set_xlim(self.resultados_simulacion[-2][0], self.resultados_simulacion[-2][-1])
            
            acum_espesor = 0
            for i in self.tab_datos_entrada.parametros.ESPESORES:
                acum_espesor += i
                self.tab_resultados_termicos.ax_flujo_calor.plot([(acum_espesor)*1000,(acum_espesor)*1000],[self.tab_resultados_termicos.min_flujo,  self.tab_resultados_termicos.max_flujo], linestyle="dashed", lw=1, color="C3")
            self.tab_resultados_termicos.slider_simulacion.configure(state=customtkinter.NORMAL)
            self.tab_resultados_termicos.canvas.draw_idle()
            Alerta(self.master.master, text="Simulación númerica exitosa", time=1000, type_of_alert=0)
        except:
            Alerta(self.master.master, text="Algo salio mal, revisar los parámetros de entrada y reintentar la simulación", time=1000, type_of_alert=-1)

    def solucion_numerica_multiple(self):
        try:
            self.resultados_simulacion_multiple = []
            # OBTENER RESULTADOS DE SIMULACIÓN MULTIPLE
            for muro in self.tab_comparacion_multiple.muros:
                self.resultados_simulacion_multiple.append(solver(espesores=muro[2].ESPESORES, rho_cp=muro[2].PRODUCTO_RHO_CP, 
                            conductividades=muro[2].CONDUCTIVIDADES, difusividades=muro[2].DIFUSIVIDADES,
                            h_i=muro[2].COEFICIENTE_CONVECCION_INTERIOR, h_e=muro[2].COEFICIENTE_CONVECCION_EXTERIOR,
                            T_e_max=muro[2].TEMPERATURA_EXTERNA_MAXIMA, rango_temp=muro[2].RANGO_TEMPERATURA,
                            T_i=muro[2].TEMPERATURA_INTERIOR, zona_horaria=muro[2].ZONA_HORARIA,
                            fecha=muro[2].FECHA, latitud=muro[2].LATITUD, 
                            orientacion=muro[2].ORIENTACION, inclinacion=muro[2].INCLINACION))
            # GANANICA Y CARGA TÉRMICA
            for line in self.tab_comparacion_multiple.ax_ganancia_termica.lines[:]:
                line.remove()
            for line in self.tab_comparacion_multiple.ax_carga_termica.lines[:]:
                line.remove()

            for i in range(len(self.tab_comparacion_multiple.muros)):
                # GRÁFICO DE GANANCIA
                self.tab_comparacion_multiple.ax_ganancia_termica.plot(self.resultados_simulacion_multiple[i][-1], self.resultados_simulacion_multiple[i][-5], color=f"C{i}", label=self.tab_comparacion_multiple.muros[i][1])
                self.tab_comparacion_multiple.ax_ganancia_termica.legend()
                
                # GRÁFICO DE CARGA
                self.tab_comparacion_multiple.ax_carga_termica.plot(self.resultados_simulacion_multiple[i][-1], self.resultados_simulacion_multiple[i][-4], color=f"C{i}", linestyle="--", label=self.tab_comparacion_multiple.muros[i][1])
                self.tab_comparacion_multiple.ax_carga_termica.legend()

            self.tab_comparacion_multiple.canvas.draw_idle()
            Alerta(self.master.master, text="Simulación númerica exitosa", time=1000, type_of_alert=0)
        except:
            Alerta(self.master.master, text="Algo salio mal, revisar los parámetros de entrada y reintentar la simulación", time=1000, type_of_alert=-1)
# APLICACIÓN
class App(customtkinter.CTk):
    def __init__(self):
        super().__init__()
        
        self.grid_rowconfigure(0, weight=1)  # configure grid system
        self.grid_columnconfigure(0, weight=1)

        self.tab_view = MyTabView(master=self)
        self.tab_view.configure(corner_radius=15, segmented_button_selected_color="#083a5e", segmented_button_selected_hover_color="#022238")
        self.tab_view.grid(row=0, column=0, padx=20, pady=10, sticky="nsew")

customtkinter.set_appearance_mode("dark")
app = App()
=======
import tkinter
import customtkinter
import json
import matplotlib.pyplot as plt
from matplotlib.backend_bases import key_press_handler
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.figure import Figure
from matplotlib.patches import Polygon, FancyArrowPatch
from matplotlib.gridspec import GridSpec
import numpy as np
import datetime
from fipy import Variable, FaceVariable, CellVariable, Grid1D, ImplicitSourceTerm, TransientTerm, DiffusionTerm, Viewer
from fipy.tools import numerix
from scipy.interpolate import CubicSpline

plt.switch_backend('agg')
DIVISION_HORA = 2

# VALIDACIONES DE PARÁMETROS DE ENTRADA    
def validacion_fecha(var):
    try:
        dia = var.timetuple().tm_yday
        mes = var.timetuple().tm_mon
        return True
    except:
        return False
def validacion_zona_horaria(var):
    try:
        if -13 < int(var) < 13:
            return True
        return False
    except:
        return False
def validacion_latitud(var):
    try:
        if -90 <= float(var) <= 90:
            return True
        return False
    except:
        return False
def validacion_orientacion(var):
    try:
        if 0 <= float(var) <= 360:
            return True
        return False
    except:
        return False
def validacion_inclinacion(var):
    try:
        if 0 <= float(var) <= 90:
            return True
        return False
    except:
        return False
def validacion_propiedades_termicas(var):
    try:
        if float(var) > 0:
            return True
        return False
    except:
        return False
def validacion_temperaturas(var):
    try:
        if float(var):
            return True
        return False
    except:
        return False
def validacion_capas(var):
    if var.isdigit():
        return True
    return False
def validacion_propiedades(var):
    if var.isdigit() or var == ".":
        return True
    return False

def solver(espesores, rho_cp, conductividades, difusividades, h_i, h_e, T_e_max, rango_temp, T_i, zona_horaria, fecha, latitud, orientacion, inclinacion):
    """
    Esta función se dedica a la solución numérica del problema
    """

    # MURO
    """
    Se define el tamaño total y el número de puntos a lo largo de este (función de las capas del muro)
    """
    L = sum(espesores) # LARGO DEL MURO
    N_PUNTOS = 1000*int(L/min(espesores)) # NÚMERO DE PUNTOS
    DX = L / N_PUNTOS # DISTANCIA ENTRE PUNTOS
    mesh = Grid1D(nx=N_PUNTOS, dx=DX) # MALLADO DEL MURO

    # TIEMPO
    """
    Se define el intervalo de tiempo y los pasos (1 hora y 24 pasos para considerar un día completo)
    """
    N_DIAS = 5
    DT = 3600/DIVISION_HORA
    PASOS = int(24*3600/DT)*N_DIAS

    # VARIABLES PARA EL MODELO
    """
    Las variables que tendremos en cuenta son el tiempo, la temperatura y la difusividad termica
    """
    tiempo = Variable() # VARIABLE TIEMPO
    temperatura = CellVariable(name="Temperatura", mesh=mesh) # VARIABLE TEMPERATURA
    D = FaceVariable(mesh=mesh) # DIFUSIVIDADES PARA CADA CARA
    K_cell = CellVariable(mesh=mesh) # CONDUCTIVIDADES PARA CADA NODO
    X = mesh.faceCenters[0] # PUNTOS DE LAS CARAS
    RHO_C = CellVariable(mesh=mesh) # PRODUCTO DENSIDAD Y CALOR ESPECIFICO PARA CADA NODO

    for j in range(N_PUNTOS ): 
        for i in range(len(espesores)):
            if X[j] >= sum(espesores[0:i]):
                D[j] = difusividades[i]
                K_cell[j] = conductividades[i]
                RHO_C[j] = conductividades[i]/difusividades[i]

    K = FaceVariable(mesh=mesh) # CONDUCTIVIDADES PARA CADA CARA
    K.setValue(K_cell.harmonicFaceValue)

    # CONDICIONES DE BORDE
    """
    Se tendrán en cuenta efectos convectivos y de radiación, temperatura externa variable y temperatura interior fija
    """
    # TEMPERATURA EXTERNA VARIABLE
    RANGO_TEMPERATURA = rango_temp
    TEMPERTURA_EXTERIOR_MAXIMA = T_e_max
    t_aprox = numerix.linspace(0,24*3600,25)
    PORCENTAJES = numerix.array([0.82, 0.87,0.92,0.96,0.99,1,0.98,0.93,0.84,0.71,0.56,0.39,0.23,0.11,0.03,0,0.03,0.10,0.21,0.34,0.47,0.58,0.68,0.76,0.82])
    cs_aprox = CubicSpline(t_aprox, PORCENTAJES)
    porcentaje = Variable()
    temperatura_exterior = TEMPERTURA_EXTERIOR_MAXIMA - RANGO_TEMPERATURA*porcentaje

    # RADIACION SOLAR
    ORIENTACION = orientacion
    LATITUD = latitud
    INCLINACION = inclinacion
    n_dia = fecha.timetuple().tm_yday
    n_mes = fecha.timetuple().tm_mon
    N_dia = numerix.deg2rad((n_dia-1)*(360/365)) 
    angulo_h = numerix.deg2rad((tiempo % (24*3600) -12*3600)/(24*3600)*360) # ÁNGULO HORARIO (UTC)
    angulo_delta = numerix.deg2rad(0.3963723-22.9132745*numerix.cos(N_dia) + 4.0254304*numerix.sin(N_dia) - 0.3872050*numerix.cos(2*N_dia) + 0.05196728*numerix.sin(2*N_dia) - 0.1545267*numerix.cos(3*N_dia) + 0.08479777*numerix.sin(3*N_dia)) # ÁNGULO DE DECLINACIÓN
    CTE_SOLAR_A = [1202, 1187, 1164, 1130, 1106, 1092, 1093, 1107, 1136, 1166, 1190, 1204]
    CTE_SOLAR_B = [0.141, 0.142, 0.149, 0.164, 0.177, 0.185, 0.186, 0.182, 0.165, 0.152, 0.142, 0.141]
    CTE_SOLAR_C = [0.103, 0.104, 0.109, 0.120, 0.130, 0.137, 0.138, 0.134, 0.121, 0.111, 0.106, 0.103]
    G_T = Variable() # RADIACION TOTAL

    # CONVECCIÓN EXTERIOR + RADIACIÓN EXTERIOR
    mask_e = mesh.facesLeft
    MA = numerix.MA
    tmp = MA.repeat(mesh._faceCenters[..., numerix.NewAxis,:], 2, 1)
    cellToFaceDistanceVectors = tmp - numerix.take(mesh._cellCenters, mesh.faceCellIDs, axis=1)
    tmp = numerix.take(mesh._cellCenters, mesh.faceCellIDs, axis=1)
    tmp = tmp[..., 1,:] - tmp[..., 0,:]
    cellDistanceVectors = MA.filled(MA.where(MA.getmaskarray(tmp), cellToFaceDistanceVectors[:, 0], tmp))
    dPf = FaceVariable(mesh=mesh, value=mesh._faceToCellDistanceRatio * cellDistanceVectors)
    n = mesh.faceNormals
    a_e = h_e*n
    b_e = conductividades[0]
    g_e = h_e*temperatura_exterior + G_T
    coef_robin_e = (mask_e * conductividades[0] * n / (-dPf.dot(a_e) + b_e))

    # TEMPERATURA INTERNA FIJA
    temperatura_interior = T_i # TEMPERATURA INTERNA FIJA

    # CONVECCION INTERIOR
    mask_i = mesh.facesRight
    MA = numerix.MA
    tmp = MA.repeat(mesh._faceCenters[..., numerix.NewAxis,:], 2, 1)
    cellToFaceDistanceVectors = tmp - numerix.take(mesh._cellCenters, mesh.faceCellIDs, axis=1)
    tmp = numerix.take(mesh._cellCenters, mesh.faceCellIDs, axis=1)
    tmp = tmp[..., 1,:] - tmp[..., 0,:]
    cellDistanceVectors = MA.filled(MA.where(MA.getmaskarray(tmp), cellToFaceDistanceVectors[:, 0], tmp))
    dPf = FaceVariable(mesh=mesh, value=mesh._faceToCellDistanceRatio * cellDistanceVectors)
    n = mesh.faceNormals
    a_i = h_i*n
    b_i = conductividades[-1]
    g_i = h_i*temperatura_interior
    coef_robin_i = (mask_i * conductividades[-1] * n / (-dPf.dot(a_i) + b_i))

    # SOLUCIÓN
    """
    La solución considera el término difusivo, el transitorio y las condiciones de borde
    """
    t = numerix.linspace(0, 24*DIVISION_HORA, 24*DIVISION_HORA + 1)
    K.setValue(0., where= mask_e | mask_i)
    eq = (TransientTerm(coeff=RHO_C) == DiffusionTerm(coeff=K) + (coef_robin_e * g_e).divergence - ImplicitSourceTerm(coeff=(coef_robin_e * h_e).divergence) 
        + (coef_robin_i * g_i).divergence - ImplicitSourceTerm(coeff=(coef_robin_i * h_i).divergence)) 

    # VARIABLAS A DEVOLVER
    temperaturas_pasos = []
    temperatura_exterior_pasos = []
    temperaturas_superficie_externa = []
    temperaturas_superficie_interna = []
    temperaturas_superficie_intermedia = []
    flujo_conduccion_entrada = []
    flujo_conduccion_salida = []
    flujo_entrada = []
    flujo_salida = []
    flujo = []

    # LOOP DE TIEMPO
    for paso in range(PASOS + 1):
        tiempo.setValue(paso * DT)  
        porcentaje.setValue(cs_aprox(paso*DT % (24*3600)))
        angulo_beta = numerix.arcsin(numerix.clip((numerix.cos(angulo_delta)*numerix.cos(LATITUD)*numerix.cos(angulo_h)+numerix.sin(angulo_delta)*numerix.sin(LATITUD)), -1, 1)) # ÁNGULO DE ALTITUD SOLAR
        if angulo_h > 0:
            angulo_phi = -numerix.arccos(numerix.clip((numerix.sin(angulo_delta)*numerix.cos(LATITUD)-numerix.sin(LATITUD)*numerix.cos(angulo_delta)*numerix.cos(angulo_h))/(numerix.cos(angulo_beta)), -1, 1)) + 2*numerix.pi
        else:
            angulo_phi = numerix.arccos(numerix.clip((numerix.sin(angulo_delta)*numerix.cos(LATITUD)-numerix.sin(LATITUD)*numerix.cos(angulo_delta)*numerix.cos(angulo_h))/(numerix.cos(angulo_beta)), -1, 1)) # ÁNGULO SOLAR AZIMUTUAL
        angulo_gamma = numerix.absolute(angulo_phi - ORIENTACION) # ÁNGULO SOLAR AZIMUTUAL RESPECTO A LA SUPERFICIE
        angulo_theta = numerix.arccos(numerix.clip(numerix.cos(angulo_beta)*numerix.cos(angulo_gamma)*numerix.sin(INCLINACION) + numerix.sin(angulo_beta)*numerix.cos(INCLINACION), -1, 1)) # ÁNGULO DE INCIDENCIA
        if numerix.sin(angulo_beta) <= .01:
            G_ND = 0
            G_T.setValue(0)
        else:
            G_ND = CTE_SOLAR_A[n_mes-1]/(numerix.exp(CTE_SOLAR_B[n_mes-1]/numerix.sin(angulo_beta)))
            if numerix.rad2deg(INCLINACION) == 90:
                if numerix.cos(angulo_theta) <= - 0.2:
                    G_T.setValue(G_ND*(max(numerix.cos(angulo_theta), 0) + 0.45*CTE_SOLAR_C[n_mes-1]))
                else:
                    G_T.setValue(G_ND*(max(numerix.cos(angulo_theta), 0) + (0.55+0.437*numerix.cos(angulo_theta)+0.313*numerix.cos(angulo_theta)**2)*CTE_SOLAR_C[n_mes-1]))
            else:
                 G_T.setValue(G_ND*(max(numerix.cos(angulo_theta), 0) + (1+numerix.cos(INCLINACION))/2*CTE_SOLAR_C[n_mes-1]))

        eq.solve(var=temperatura, dt=DT)
        temperaturas_pasos.append(temperatura.faceValue.copy())
        temperatura_exterior_pasos.append(temperatura_exterior.copy())
        temperaturas_superficie_externa.append(temperaturas_pasos[paso][0])
        temperaturas_superficie_interna.append(temperaturas_pasos[paso][-1])
        temperaturas_superficie_intermedia.append(temperaturas_pasos[paso][int(len(X)/2)])
        flujo_conduccion_entrada.append(-conductividades[0]*temperatura.faceGrad[0][1].copy())
        flujo_conduccion_salida.append(-conductividades[-1]*temperatura.faceGrad[0][-2].copy())
        flujo_entrada.append(h_e*(temperatura_exterior.copy() - temperaturas_superficie_externa[paso]) + G_T.copy())
        flujo_salida.append(h_i*(-temperatura_interior + temperaturas_superficie_interna[paso]))
        flujo_paso = FaceVariable(mesh=mesh,value=-K*temperatura.faceGrad[0].copy())
        flujo_paso[0] = flujo_paso[1]
        flujo_paso[-1] = flujo_paso[-2]
        flujo.append(flujo_paso) # CORRECCIÓN DEL FLUJO

    # SE RECORTAN LOS VALORES
    temperatura_exterior_pasos = temperatura_exterior_pasos[-24*DIVISION_HORA-1:]
    temperaturas_pasos = temperaturas_pasos[-24*DIVISION_HORA-1:]
    temperaturas_superficie_externa = temperaturas_superficie_externa[-24*DIVISION_HORA-1:]
    temperaturas_superficie_interna = temperaturas_superficie_interna[-24*DIVISION_HORA-1:]
    temperaturas_superficie_intermedia = temperaturas_superficie_intermedia[-24*DIVISION_HORA-1:]
    flujo_conduccion_entrada = flujo_conduccion_entrada[-24*DIVISION_HORA-1:]
    flujo_conduccion_salida = flujo_conduccion_salida[-24*DIVISION_HORA-1:]
    flujo_entrada = flujo_entrada[-24*DIVISION_HORA-1:]
    flujo_salida = flujo_salida[-24*DIVISION_HORA-1:]
    flujo = flujo[-24*DIVISION_HORA-1:]

    
    return [temperaturas_pasos, temperaturas_superficie_externa, temperaturas_superficie_intermedia, temperaturas_superficie_interna, flujo_conduccion_entrada, 
            flujo_conduccion_salida, flujo, X*1000, t/DIVISION_HORA]


# ALERTA
class Alerta(customtkinter.CTkLabel):
    def __init__(self, master, text, type_of_alert,  time):
        super().__init__(master)
        if type_of_alert == 0:
            color = "#007a60"
        else:
            color = "#6e0505"
        self.configure(text=text, text_color="white", fg_color=color)
        self.place(relwidth=1, relx=0, rely=0, anchor="nw", relheight=0.05)
        self.after(time, lambda: self.destroy())

# DATOS DE ENTRADA
class Parametros_entrada:
    def __init__(self):
        
        # CONDICIONES GEOGRÁFICAS
        self.FECHA = 0 # DIA Y MES
        self.ZONA_HORARIA = 0 
        self.LATITUD = 0 
        self.ORIENTACION = 0 # ORIENTACIÓN DEL MURO RESPECTO AL NORTE
        self.INCLINACION = 0 # SUPERFICIE VERTICAL, INCLINADA U HORIZONTAL
        # CONDICIONES INTERNAS
        self.TEMPERATURA_INTERIOR = 0

        # CONDICIONES EXTERNAS
        self.TEMPERATURA_EXTERNA_MAXIMA = 0
        self.RANGO_TEMPERATURA = 0

        # DEFINICIÓN MURO
        self.N_CAPAS = 0
        self.ESPESORES = []
        self.DENSIDADES = []
        self.CALORES_ESPECIFICOS = []
        self.CONDUCTIVIDADES = []
        self.DIFUSIVIDADES = []
        self.PRODUCTO_RHO_CP = []

        # CONVECCIÓN INTERIOR
        self.COEFICIENTE_CONVECCION_INTERIOR = 0

        # CONVECCIÓN EXTERIOR
        self.COEFICIENTE_CONVECCION_EXTERIOR = 0

# PESTAÑA DATOS DE ENTRADA
class Tab_datos_entrada(customtkinter.CTkFrame):
    def __init__(self, master, controller):
        super().__init__(master)
        self.controller = controller
        self.datos_guardados = False
        self.parametros = Parametros_entrada()
        master.grid_columnconfigure(0, weight = 1)
        self.grid_columnconfigure((0,1,2,3,4,5), weight=1)
        
        # FECHA
        self.label_fecha = customtkinter.CTkLabel(self, text="FECHA: ")
        self.label_fecha.grid(row=0, column=0, padx=20, pady=0, sticky="w")
        self.fecha = customtkinter.CTkEntry(self, placeholder_text="DD/MM")
        self.fecha.grid(row=0, column=1, padx=20, pady=0, sticky="ew")

        # ZONA HORARIA
        self.label_zona_horaria = customtkinter.CTkLabel(self, text="ZONA HORARIA: (RESPECTO A UTC)")
        self.label_zona_horaria.grid(row=1, column=0, padx=20, pady=0, sticky="w")
        self.zona_horaria = customtkinter.CTkEntry(self, placeholder_text="-12 a 12")
        self.zona_horaria.grid(row=1, column=1, padx=20, pady=0, sticky="ew")

        # LATITUD
        self.label_latitud = customtkinter.CTkLabel(self, text="LATITUD:")
        self.label_latitud.grid(row=2, column=0, padx=20, pady=0, sticky="w")
        self.latitud = customtkinter.CTkEntry(self, placeholder_text="- 90° a 90°")
        self.latitud.grid(row=2, column=1, padx=20, pady=0, sticky="ew")

        # ORIENTACION
        self.label_orientacion = customtkinter.CTkLabel(self, text="ORIENTACIÓN MURO: ")
        self.label_orientacion.grid(row=3, column=0, padx=20, pady=0, sticky="w")
        self.orientacion = customtkinter.CTkEntry(self, placeholder_text="0° a 360°")
        self.orientacion.grid(row=3, column=1, padx=20, pady=0, sticky="ew")

        # INCLINACION
        self.label_inclinacion = customtkinter.CTkLabel(self, text="INCLINACIÓN MURO: ")
        self.label_inclinacion.grid(row=4, column=0, padx=20, pady=0, sticky="w")
        self.inclinacion = customtkinter.CTkEntry(self, placeholder_text="0° a 90°")
        self.inclinacion.grid(row=4, column=1, padx=20, pady=0, sticky="ew")

        # MAPA
        fig = Figure(facecolor='#2b2b2b')
        ax = fig.add_subplot()
        ax.axis('off')
        fig.tight_layout()
        fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
        img = plt.imread("imagen.png")
        img_altura = 1488
        img_ancho = 2976
        ax.imshow(img, extent=[0,img_ancho,0,img_altura])
        ax.set_xlim(0, img_ancho)
        ax.set_ylim(0, img_altura)
        canvas = FigureCanvasTkAgg(fig, master=self) 
        canvas.draw()
        canvas.get_tk_widget().grid(row=0, column=2, rowspan=6, columnspan=4, sticky="nsew")
        poligono_zona_horaria = None
        vector = None
        muro = None
        # RECARGAR EL MAPA CON LOS DATOS CORRECTOS
        def recargar_mapa(poligono_zona_horaria, vector, muro):
            try:
                for patch in ax.patches[:]:  # Use a copy of the list to avoid modifying it while iterating
                    patch.remove()
                for line in ax.lines[:]:
                    line.remove()
                fecha = datetime.datetime.strptime(self.fecha.get() + "/25", "%d/%m/%y")
                latitud = float(self.latitud.get())
                zona_horaria = int(self.zona_horaria.get())
                orientacion = float(self.orientacion.get())
                inclinacion = float(self.inclinacion.get())
                if not validacion_fecha(fecha):
                    raise ValueError("Fecha no valida")
                if not validacion_latitud(latitud):
                    raise ValueError("Latitud no valida")
                if not validacion_zona_horaria(zona_horaria):
                    raise ValueError("Zona horaria no valida")
                if not validacion_orientacion(orientacion): 
                    raise ValueError("Orientacion no valida")
                if not validacion_inclinacion(inclinacion):
                      raise ValueError("Inclinacion no valida")
                
                if zona_horaria == -12:
                    inicio_longitud = zona_horaria*15
                    fin_longitud = zona_horaria*15 + 7.5
                elif zona_horaria == 12:
                    inicio_longitud = zona_horaria*15 -7.5
                    fin_longitud = zona_horaria*15
                else:
                    inicio_longitud = zona_horaria*15 - 7.5
                    fin_longitud = zona_horaria*15 + 7.5

                # POLIGONO ZONA HORARIA 
                x1, y1 = (inicio_longitud/180*img_ancho/2+img_ancho/2, -90/90*img_altura/2+img_altura/2)
                x2, y2 = (fin_longitud/180*img_ancho/2+img_ancho/2, 90/90*img_altura/2+img_altura/2)
                if poligono_zona_horaria:
                    poligono_zona_horaria.remove()
                poligono_zona_horaria = Polygon(
                    [(x1, y1), (x1, y2), (x2, y2), (x2, y1)], 
                    facecolor="#05646b", alpha=0.3, zorder=2
                )
                ax.add_patch(poligono_zona_horaria)
                
                # LATITUD
                ax.plot([0,img_ancho], [latitud/90*img_altura/2+img_altura/2,latitud/90*img_altura/2+img_altura/2], color="#8ccffc", linewidth=2)

                # MURO Y SU ORIENTACIÓN
                mid_longitud = (inicio_longitud + fin_longitud) / 2
                x0, y0 = (mid_longitud/180*img_ancho/2+img_ancho/2, latitud/90*img_altura/2+img_altura/2)
                angulo_radianes = np.radians(orientacion)
                largo_vector = 150 
                largo_muro = 50
                vector_x = x0 + largo_vector * np.sin(angulo_radianes)
                vector_y = y0 + largo_vector * np.cos(angulo_radianes)
                angulo_radianes_perpendicular = angulo_radianes + np.pi / 2 
                muro_x = x0 + largo_muro * np.sin(angulo_radianes_perpendicular)
                muro_y = y0 + largo_muro * np.cos(angulo_radianes_perpendicular)
                if vector:
                    vector.remove()
                if muro:
                    muro.remove()
                vector = FancyArrowPatch(
                    (x0, y0), (vector_x, vector_y), mutation_scale=15, color="#52bff1", lw=1, zorder=10)
                ax.add_patch(vector)
                muro = FancyArrowPatch(
                    (2*x0 - muro_x, 2*y0 - muro_y), (muro_x, muro_y), arrowstyle="-", mutation_scale=25, color="#1C5B80", lw=4, zorder=11)
                ax.add_patch(muro)
                fig.canvas.draw_idle()

            except ValueError as e:
                Alerta(self.master.master.master, text="Also salió mal, revisar los parámetros de entrada geográficos", type_of_alert=-1, time=1000)
        
        # GUARDADO DEL MAPA
        self.guardado_mapa = customtkinter.CTkButton(self, text="Mostrar en mapa", command=lambda: recargar_mapa(poligono_zona_horaria, vector, muro))
        self.guardado_mapa.grid(row=5, column=0, padx=20, pady=5, sticky="ew", columnspan=2)

        # TEMPERATURA EXTERIOR
        self.label_temperatura_exterior = customtkinter.CTkLabel(self, text="TEMPERATURA EXTERIOR MÁXIMA: ")
        self.label_temperatura_exterior.grid(row=6, column=0, padx=20, pady=10, sticky="w")
        self.temperatura_exterior_maxima = customtkinter.CTkEntry(self, placeholder_text="°C")
        self.temperatura_exterior_maxima.grid(row=6, column=1, padx=20, pady=10, sticky="ew")
        self.label_rango_temperatura = customtkinter.CTkLabel(self, text="RANGO DE TEMPERATURA EXTERIOR: ")
        self.label_rango_temperatura.grid(row=6, column=2, padx=20, pady=10, sticky="w")
        self.rango_temperatura = customtkinter.CTkEntry(self, placeholder_text="°C")
        self.rango_temperatura.grid(row=6, column=3, padx=20, pady=10, sticky="ew")

        # CONVECCION EXTERIOR
        self.label_conveccion_exterior = customtkinter.CTkLabel(self, text="COEFICIENTE CONVECCION EXTERIOR: ")
        self.label_conveccion_exterior.grid(row=6, column=4, padx=20, pady=10, sticky="w")
        self.coeficiente_conveccion_exterior = customtkinter.CTkEntry(self, placeholder_text="W/m²K")
        self.coeficiente_conveccion_exterior.grid(row=6, column=5, padx=20, pady=10, sticky="ew")

        # TEMPERATURA INTERIOR
        self.label_temperatura_interior = customtkinter.CTkLabel(self, text="TEMPERATURA INTERIOR: ")
        self.label_temperatura_interior.grid(row=7, column=0, padx=20, pady=10, sticky="w")
        self.temperatura_interior = customtkinter.CTkEntry(self, placeholder_text="°C")
        self.temperatura_interior.grid(row=7, column=1, padx=20, pady=10, sticky="ew")

        # CONVECCION INTERIOR
        self.label_conveccion_interior = customtkinter.CTkLabel(self, text="COEFICIENTE CONVECCIÓN INTERIOR: ")
        self.label_conveccion_interior.grid(row=7, column=2, padx=20, pady=10, sticky="w")
        self.coeficiente_conveccion_interior = customtkinter.CTkEntry(self, placeholder_text="W/m²K")
        self.coeficiente_conveccion_interior.grid(row=7, column=3, padx=20, pady=10, sticky="ew")

        # NÚMERO DE CAPAS DEL MURO
        self.label_n_capas = customtkinter.CTkLabel(self, text="NÚMERO DE CAPAS: ")
        self.label_n_capas.grid(row=7, column=4, padx=20, pady=10, sticky="w")
        self.n_capas = customtkinter.CTkEntry(self, placeholder_text="", validate="key", validatecommand=(self.register(validacion_capas), "%S"))
        self.n_capas.grid(row=7, column=5, padx=20, pady=10, sticky="ew")
        self.propiedades_muro = ["N° DE CAPA", "ESPESOR [mm]", "DENSIDAD [kg/m³]", "CALOR ESPECÍFICO [J/kgK]", "CONDUCTIVIDAD TÉRMICA [W/mK]"]
        self.labels_temporales = []
        self.labels_capas = []
        self.espesores_capas = []
        self.densidades_capas = []
        self.calores_especificos_capas = []
        self.conductividades_termicas_capas = []

        self.boton_n_capas = customtkinter.CTkButton(self, text="Seleccionar capas", command=self.cargar_capas)
        self.boton_n_capas.grid(row=8, column=0, padx=20, pady=10, sticky="ew", columnspan=6)

        self.guardar_todo = customtkinter.CTkButton(self, text="Guardar todo", fg_color="#23a15c", hover_color="#136337", command=self.carga_de_datos)

    # TABLA CON DATOS DEL MURO
    def cargar_capas(self):
        if self.n_capas.get() == "":
            return True
        else:
            n_capas = max(1, min(int(self.n_capas.get()), 5))
            if len(self.labels_temporales) == 0:
                for i in range(len(self.propiedades_muro)):
                    tmp = customtkinter.CTkLabel(self, text=self.propiedades_muro[i])
                    tmp.grid(row = 9 , column=i, padx=20, pady=10)
                    self.labels_temporales.append(tmp)
            else:
                for i in range(len(self.labels_capas)):
                    self.labels_capas[i].destroy()
                    self.espesores_capas[i].destroy()
                    self.densidades_capas[i].destroy()
                    self.calores_especificos_capas[i].destroy()
                    self.conductividades_termicas_capas[i].destroy()
                self.labels_capas = []
                self.espesores_capas = []
                self.densidades_capas = []
                self.calores_especificos_capas = []
                self.conductividades_termicas_capas = []
                self.guardar_todo.grid_forget()
            for i in range(n_capas):
                    label_capa = customtkinter.CTkLabel(self, text=i + 1)
                    label_capa.grid(row = 10 + i, column= 0, padx=20, pady=10)
                    self.labels_capas.append(label_capa)
                    label_espesor = customtkinter.CTkEntry(self, validate="key", validatecommand=(self.register(validacion_propiedades), "%S"))
                    label_espesor.grid(row = 10 + i, column = 1, padx=20, pady=10)
                    self.espesores_capas.append(label_espesor)
                    label_densidad = customtkinter.CTkEntry(self, validate="key", validatecommand=(self.register(validacion_propiedades), "%S"))
                    label_densidad.grid(row = 10 + i, column = 2, padx=20, pady=10)
                    self.densidades_capas.append(label_densidad)
                    label_calor_especifico = customtkinter.CTkEntry(self, validate="key", validatecommand=(self.register(validacion_propiedades), "%S"))
                    label_calor_especifico.grid(row = 10 + i, column = 3, padx=20, pady=10)
                    self.calores_especificos_capas.append(label_calor_especifico)
                    label_conductividad = customtkinter.CTkEntry(self, validate="key", validatecommand=(self.register(validacion_propiedades), "%S"))
                    label_conductividad.grid(row = 10 + i, column = 4, padx=20, pady=10)
                    self.conductividades_termicas_capas.append(label_conductividad)
            self.guardar_todo.grid(row=10, column=5, padx=20, pady=10, rowspan=n_capas, sticky="nsew")

    # CARGAR LOS DATOS Y VALIDAR LAS ENTRADAS
    def carga_de_datos(self):
        try:
            self.parametros.FECHA = datetime.datetime.strptime(self.fecha.get() + "/25", "%d/%m/%y")
            self.parametros.LATITUD = numerix.deg2rad(float(self.latitud.get()))
            self.parametros.ZONA_HORARIA = numerix.deg2rad(int(self.zona_horaria.get()))
            self.parametros.ORIENTACION = numerix.deg2rad(float(self.orientacion.get()))
            self.parametros.INCLINACION = numerix.deg2rad(float(self.inclinacion.get()))
            if (not validacion_fecha(self.parametros.FECHA) or not validacion_latitud(self.parametros.LATITUD) or not validacion_zona_horaria(self.parametros.ZONA_HORARIA) or
                not validacion_orientacion(self.parametros.ORIENTACION) or not validacion_inclinacion(self.parametros.INCLINACION)):
                raise ValueError("Datos geográficos inválidos")
            self.parametros.COEFICIENTE_CONVECCION_INTERIOR = float(self.coeficiente_conveccion_interior.get())
            self.parametros.COEFICIENTE_CONVECCION_EXTERIOR = float(self.coeficiente_conveccion_exterior.get())
            self.parametros.TEMPERATURA_INTERIOR = float(self.temperatura_interior.get())
            self.parametros.TEMPERATURA_EXTERNA_MAXIMA = float(self.temperatura_exterior_maxima.get())
            self.parametros.RANGO_TEMPERATURA = float(self.rango_temperatura.get())
            self.parametros.N_CAPAS = max(1, min(int(self.n_capas.get()), 5))
            if (not validacion_propiedades_termicas(self.parametros.COEFICIENTE_CONVECCION_INTERIOR) or not validacion_propiedades_termicas(self.parametros.COEFICIENTE_CONVECCION_EXTERIOR)
            or not validacion_temperaturas(self.parametros.TEMPERATURA_INTERIOR) or not validacion_temperaturas(self.parametros.TEMPERATURA_EXTERNA_MAXIMA) or
            not validacion_temperaturas(self.parametros.RANGO_TEMPERATURA)):
                raise ValueError("Datos térmicos inválidos")
            if len(self.parametros.ESPESORES) != 0:
                self.parametros.ESPESORES = []
                self.parametros.DENSIDADES = []
                self.parametros.CALORES_ESPECIFICOS = []
                self.parametros.CONDUCTIVIDADES = []
                self.parametros.DIFUSIVIDADES = []
                self.parametros.PRODUCTO_RHO_CP = []
            for i in range(self.parametros.N_CAPAS):
                self.parametros.ESPESORES.append(float(self.espesores_capas[i].get())/1000)
                self.parametros.DENSIDADES.append(float(self.densidades_capas[i].get()))
                self.parametros.CALORES_ESPECIFICOS.append(float(self.calores_especificos_capas[i].get()))
                self.parametros.CONDUCTIVIDADES.append(float(self.conductividades_termicas_capas[i].get()))
                self.parametros.DIFUSIVIDADES.append(self.parametros.CONDUCTIVIDADES[i]/(self.parametros.CALORES_ESPECIFICOS[i]*self.parametros.DENSIDADES[i]))
                self.parametros.PRODUCTO_RHO_CP.append(self.parametros.CALORES_ESPECIFICOS[i]*self.parametros.DENSIDADES[i])
                if (not validacion_propiedades_termicas(self.parametros.ESPESORES[i]) or not validacion_propiedades_termicas(self.parametros.DENSIDADES[i])
                or not validacion_propiedades_termicas(self.parametros.CALORES_ESPECIFICOS[i]) or not validacion_propiedades_termicas(self.parametros.CONDUCTIVIDADES[i])):
                    raise ValueError("Datos del muro inválidos")
            self.datos_guardados = True
            self.controller.tab_perfiles_temperatura.boton_simulacion.configure(state=customtkinter.NORMAL)
            self.controller.tab_resultados_termicos.boton_simulacion.configure(state=customtkinter.NORMAL)
            Alerta(self.master.master.master, text="Datos guardados correctamente", type_of_alert=0, time=1000)
        except ValueError as e:
            if str(e) == "Datos geográficos inválidos" or str(e) == "Datos térmicos inválidos":
                Alerta(self.master.master.master, text=str(e), type_of_alert=-1, time=1000)
            else:
                Alerta(self.master.master.master, text="Algo salio mal, revisar los parámetros de entrada", type_of_alert=-1, time=1000)
        except:
            Alerta(self.master.master.master, text="Algo salio mal, revisar los parámetros de entrada", type_of_alert = -1, time=1000)
    
            
# PESTAÑA PERFILES DE TEMPERATURA
class Tab_perfiles_de_temperaturas(customtkinter.CTkFrame):
    def __init__(self, master, controller):
        super().__init__(master)
        self.controller = controller
        master.grid_columnconfigure(0, weight = 1)
        master.grid_rowconfigure(0, weight = 1)
        self.grid_columnconfigure((0,1,2,3,4,5), weight=1)
        self.grid_rowconfigure((0,1,2,3,4,5), weight=1)

        # GRÁFICOS
        plt.style.use("dark_background")
        self.fig = Figure(facecolor='#2b2b2b')
        gs = GridSpec(2, 1, wspace=0.2, hspace=0.4, left=0.1, right=0.9, top=0.95, bottom=0.06)

        # GRÁFICO TEMPERATURA DEL MURO
        self.ax_temp_muro = self.fig.add_subplot(gs[0, 0]) 
        self.ax_temp_muro.set_title("Perfil de temperatura en el muro")
        self.ax_temp_muro.set_xlabel(r"$e$ [mm]")
        self.ax_temp_muro.set_ylabel(r"$T$ [°C]")

        # GRÁFICO TEMPERATURAS EN SUPERFICIES
        self.ax_temp_superficie = self.fig.add_subplot(gs[1, 0])
        self.ax_temp_superficie.set_xticks([0,2,4,6,8,10,12,14,16,18,20,22,24])
        self.ax_temp_superficie.set_title("Temperatura en superficie exterior, interior e intermedia")
        self.ax_temp_superficie.set_xlim(0,24)
        self.ax_temp_superficie.set_xlabel(r"Hora [h]")
        self.ax_temp_superficie.set_ylabel(r"$T$ [°C]")
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=self) 
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=0, column=0, rowspan=5, columnspan=6, sticky="nsew")
        
        # CONTROLES
        self.label_slider = customtkinter.CTkLabel(self, text="HORA A ANALIZAR: 00:00")
        self.label_slider.grid(row=5,column=0, sticky="e", padx=20, pady=10)
        self.slider_simulacion = customtkinter.CTkSlider(self, from_=0, to=24, number_of_steps=24, state=customtkinter.DISABLED, command=self.seleccionar_hora)
        self.slider_simulacion.set(0)
        self.slider_simulacion.grid(row=5,column=1, columnspan=2, sticky="ew", padx=20, pady=10)
        self.boton_simulacion = customtkinter.CTkButton(self, text="Correr simulación", command=self.controller.solucion_numerica, state=customtkinter.DISABLED)
        self.boton_simulacion.grid(row=5, column=4, sticky="ew", padx=20, pady=10)

        self.min_temp = 0
        self.max_temp = 0

    def seleccionar_hora(self, value):
        value = int(value)
        for line in self.ax_temp_muro.lines[:]:
            line.remove()
        self.ax_temp_muro.plot(self.controller.resultados_simulacion[-2], self.controller.resultados_simulacion[0][value*DIVISION_HORA], color="C1")
        if int(value) == 0:
            value = "00"
        self.label_slider.configure(text=f"HORA A ANALIZAR: {int(value)}:00")
        acum_espesor = 0
        for i in self.controller.tab_datos_entrada.parametros.ESPESORES:
            acum_espesor += i
            self.ax_temp_muro.plot([(acum_espesor)*1000,(acum_espesor)*1000], [self.min_temp, self.max_temp], linestyle="dashed", lw=1, color="C3")
        self.canvas.draw_idle()

# PESTAÑA RESULTADOS TÉRMICOS
class Tab_resultados_termicos(customtkinter.CTkFrame):
    def __init__(self, master, controller):
        super().__init__(master)
        self.controller = controller
        master.grid_columnconfigure(0, weight = 1)
        master.grid_rowconfigure(0, weight = 1)
        self.grid_columnconfigure((0,1,2,3,4,5), weight=1)
        self.grid_rowconfigure((0,1,2,3,4,5), weight=1)
        
        # GRÁFICOS
        plt.style.use("dark_background")
        self.fig = Figure(facecolor='#2b2b2b')
        gs = GridSpec(2, 1, wspace=0.2, hspace=0.4, left=0.1, right=0.9, top=0.95, bottom=0.06)

        # GRÁFICO FLUJO DE CALOR EN EL MURO
        self.ax_flujo_calor = self.fig.add_subplot(gs[0, 0]) 
        self.ax_flujo_calor.set_title("Flujo de calor a lo largo del muro")
        self.ax_flujo_calor.set_xlabel(r"$e$ [mm]")
        self.ax_flujo_calor.set_ylabel(r"$Q$ $[W/m^2]$")

        # GRÁFICO GANANCIA y CARGA TÉRMICA
        self.ax_ganancia_carga_termica = self.fig.add_subplot(gs[1, 0])
        self.ax_ganancia_carga_termica.set_title("Ganancia y carga térmica")
        self.ax_ganancia_carga_termica.set_xlim(0,24)
        self.ax_ganancia_carga_termica.set_xticks([0,2,4,6,8,10,12,14,16,18,20,22,24])
        self.ax_ganancia_carga_termica.set_xlabel(r"Hora [h]")
        self.ax_ganancia_carga_termica.set_ylabel(r"$Q$ $[W/m^2]$")
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=self) 
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=0, column=0, rowspan=5, columnspan=6, sticky="nsew")

        # CONTROLES
        self.label_slider = customtkinter.CTkLabel(self, text="HORA A ANALIZAR: 00:00")
        self.label_slider.grid(row=5,column=0, sticky="e", padx=20, pady=10)
        self.slider_simulacion = customtkinter.CTkSlider(self, from_=0, to=24, number_of_steps=24, state=customtkinter.DISABLED, command=self.seleccionar_hora)
        self.slider_simulacion.set(0)
        self.slider_simulacion.grid(row=5,column=1, columnspan=2, sticky="ew", padx=20, pady=10)
        self.boton_simulacion = customtkinter.CTkButton(self, text="Correr simulación", command=self.controller.solucion_numerica, state=customtkinter.DISABLED)
        self.boton_simulacion.grid(row=5, column=4, sticky="ew", padx=20, pady=10)

        self.min_flujo = 0
        self.max_flujo = 0

    def seleccionar_hora(self, value):
        value = int(value)
        for line in self.ax_flujo_calor.lines[:]:
            line.remove()
        self.ax_flujo_calor.plot(self.controller.resultados_simulacion[-2], self.controller.resultados_simulacion[-3][value*DIVISION_HORA], color="C1")
        if int(value) == 0:
            value = "00"
        self.label_slider.configure(text=f"HORA A ANALIZAR: {int(value)}:00")
        acum_espesor = 0
        for i in self.controller.tab_datos_entrada.parametros.ESPESORES:
            acum_espesor += i
            self.ax_flujo_calor.plot([(acum_espesor)*1000,(acum_espesor)*1000],[self.min_flujo, self.max_flujo], linestyle="dashed", lw=1, color="C3")
        self.canvas.draw_idle()

# PESTAÑA COMPARACIÓN MÚLTIPLE
class Tab_comparacion_multiple(customtkinter.CTkFrame):
    def __init__(self, master, controller):
        super().__init__(master)
        self.controller = controller
        master.grid_columnconfigure(0, weight = 1)
        master.grid_rowconfigure(0, weight = 1)
        self.grid_columnconfigure((0,1,2,3,4,5,6), weight=1)
        self.grid_rowconfigure((0,1,2,3,4,5), weight=1)
        self.muros = []

        # BOTÓN CARGA ARCHIVO
        self.boton_carga_archivo = customtkinter.CTkButton(self, text="Cargar archivo .json", fg_color="#23a15c", hover_color="#136337", command=self.cargar_json)
        self.boton_carga_archivo.grid(row=0, column=0, padx=10, sticky="ew")

        # TEXTBOX CON DATOS DE LOS MUROS
        self.texto_json = customtkinter.CTkTextbox(self, width=550, height=300)
        self.texto_json.grid(row=1, column=0, padx=10, sticky="nsew", rowspan=4)

        # CORRER SIMULACION
        self.boton_simulacion = customtkinter.CTkButton(self, text="Correr simulación", command=self.controller.solucion_numerica_multiple, state=customtkinter.DISABLED)
        self.boton_simulacion.grid(row=5, column=0, sticky="ew", padx=10)
        
        # GRÁFICOS
        plt.style.use("dark_background")
        self.fig = Figure(facecolor='#2b2b2b')
        gs = GridSpec(2, 1, wspace=0.2, hspace=0.4, left=0.1, right=0.9, top=0.95, bottom=0.06)


        # GRÁFICO GANANCIA TÉRMICA
        self.ax_ganancia_termica = self.fig.add_subplot(gs[0, 0])
        self.ax_ganancia_termica.set_title("Ganancia térmica")
        self.ax_ganancia_termica.set_xlim(0,24)
        self.ax_ganancia_termica.set_xticks([0,2,4,6,8,10,12,14,16,18,20,22,24])
        self.ax_ganancia_termica.set_xlabel(r"Hora [h]")
        self.ax_ganancia_termica.set_ylabel(r"$Q$ $[W/m^2]$")

        # GRÁFICO CARGA TÉRMICA
        self.ax_carga_termica = self.fig.add_subplot(gs[1, 0])
        self.ax_carga_termica.set_title("Carga térmica")
        self.ax_carga_termica.set_xlim(0,24)
        self.ax_carga_termica.set_xticks([0,2,4,6,8,10,12,14,16,18,20,22,24])
        self.ax_carga_termica.set_xlabel(r"Hora [h]")
        self.ax_carga_termica.set_ylabel(r"$Q$ $[W/m^2]$")
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=self) 
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=0, column=1, rowspan=6, columnspan=5, sticky="nsew")

    # CARGA Y PROCESAMIENTO DEL ARCHIVO .JSON
    def cargar_json(self):
        archivo = tkinter.filedialog.askopenfilename(
            title="Seleccionar archivo JSON",
            filetypes=[("Archivos JSON", "*.json")]
        )
        self.muros = []
        if archivo:
            try:
                with open(archivo, "r", encoding="utf-8") as f:
                    datos = json.load(f)
                    for muro in datos:
                        id_muro = muro["id"]
                        nombre_muro = muro["nombre"]
                        parametros = Parametros_entrada()
                        parametros.FECHA = datetime.datetime.strptime(muro["fecha"] + "/25", "%d/%m/%y")
                        parametros.ZONA_HORARIA = numerix.deg2rad(muro["zona_horaria"])
                        parametros.LATITUD = numerix.deg2rad(muro["latitud"])
                        parametros.ORIENTACION = numerix.deg2rad(muro["orientacion"])
                        parametros.INCLINACION = numerix.deg2rad(muro["inclinacion"])
                        if (not validacion_fecha(parametros.FECHA) or not validacion_latitud(parametros.LATITUD) or not validacion_zona_horaria(parametros.ZONA_HORARIA) or
                            not validacion_orientacion(parametros.ORIENTACION) or not validacion_inclinacion(parametros.INCLINACION)):
                            raise ValueError("Datos geográficos inválidos")
                        parametros.COEFICIENTE_CONVECCION_INTERIOR = muro["coeficiente_conveccion_interior"]
                        parametros.COEFICIENTE_CONVECCION_EXTERIOR = muro["coeficiente_conveccion_exterior"]
                        parametros.TEMPERATURA_INTERIOR = muro["temperatura_interior"]
                        parametros.TEMPERATURA_EXTERNA_MAXIMA = muro["temperatura_exterior"]
                        parametros.RANGO_TEMPERATURA = muro["rango_temperatura_exterior"]
                        parametros.N_CAPAS = muro["n_capas"]
                        
                        if (not validacion_propiedades_termicas(parametros.COEFICIENTE_CONVECCION_INTERIOR) or not validacion_propiedades_termicas(parametros.COEFICIENTE_CONVECCION_EXTERIOR)
                        or not validacion_temperaturas(parametros.TEMPERATURA_INTERIOR) or not validacion_temperaturas(parametros.TEMPERATURA_EXTERNA_MAXIMA) or
                        not validacion_temperaturas(parametros.RANGO_TEMPERATURA)):
                            raise ValueError("Datos térmicos inválidos")
                        for capa in muro["capas"]:
                            parametros.ESPESORES.append(capa["espesor"]/1000)
                            parametros.DENSIDADES.append(capa["densidad"])
                            parametros.CALORES_ESPECIFICOS.append(capa["calor_especifico"])
                            parametros.CONDUCTIVIDADES.append(capa["conductividad"])
                            parametros.DIFUSIVIDADES.append(capa["conductividad"]/(capa["densidad"]*capa["calor_especifico"]))
                            parametros.PRODUCTO_RHO_CP.append(capa["densidad"]*capa["calor_especifico"])
                            if (not validacion_propiedades_termicas(parametros.ESPESORES[-1]) or not validacion_propiedades_termicas(parametros.DENSIDADES[-1])
                            or not validacion_propiedades_termicas(parametros.CALORES_ESPECIFICOS[-1]) or not validacion_propiedades_termicas(parametros.CONDUCTIVIDADES[-1])):
                                raise ValueError("Datos del muro inválidos")
                        self.muros.append([id_muro, nombre_muro, parametros])
                    contenido = json.dumps(datos, indent=4, ensure_ascii=False)
                    self.texto_json.delete("0.0", "end")
                    self.texto_json.insert("0.0", contenido)
                    self.boton_simulacion.configure(state=customtkinter.NORMAL)
                    Alerta(self.master.master.master, text="Archivo .json cargado correctamente", type_of_alert=0, time=1000)
            except Exception as e:
                Alerta(self.master.master.master, text="Error al cargar el archivo .json, revisar los datos ingresados", type_of_alert=-1, time=2000)

# PESTAÑAS
class MyTabView(customtkinter.CTkTabview):
    def __init__(self, master, **kwargs):
        super().__init__(master, **kwargs)
        self.resultados_simulacion = []
        self.resultados_simulacion_multiple = []

        self.grid_columnconfigure(0, weight=1)
        # PESTAÑAS
        self.add("Datos de entrada")
        self.add("Resultados temperatura")
        self.add("Resultados térmicos")
        self.add("Comparación múltiple")

        # WIDGETS DATOS DE ENTRADA
        self.tab_datos_entrada = Tab_datos_entrada(master=self.tab("Datos de entrada"), controller=self)
        self.tab_datos_entrada.grid(row=0, column=0, padx=20, pady=10, sticky="nsew")

        # WIDGETS PERFILES DE TEMPERATURA
        self.tab_perfiles_temperatura = Tab_perfiles_de_temperaturas(master=self.tab("Resultados temperatura"), controller=self)
        self.tab_perfiles_temperatura.grid(row=0, column=0, padx=20, pady=10, sticky="nsew")

        # WIDGETS RESULTADOS TÉRMICOS
        self.tab_resultados_termicos = Tab_resultados_termicos(master=self.tab("Resultados térmicos"), controller=self)
        self.tab_resultados_termicos.grid(row=0, column=0, padx=20, pady=10, sticky="nsew")

        # WIDGETS COMPARACION MULTIPLE
        self.tab_comparacion_multiple = Tab_comparacion_multiple(master=self.tab("Comparación múltiple"), controller=self)
        self.tab_comparacion_multiple.grid(row=0, column=0, padx=20, pady=10, sticky="nsew")

    def solucion_numerica(self):
        try:
            # OBTENER RESULTADOS DE SIMULACIÓN
            self.resultados_simulacion = solver(espesores=self.tab_datos_entrada.parametros.ESPESORES, rho_cp=self.tab_datos_entrada.parametros.PRODUCTO_RHO_CP, 
                        conductividades=self.tab_datos_entrada.parametros.CONDUCTIVIDADES, difusividades=self.tab_datos_entrada.parametros.DIFUSIVIDADES,
                        h_i=self.tab_datos_entrada.parametros.COEFICIENTE_CONVECCION_INTERIOR, h_e=self.tab_datos_entrada.parametros.COEFICIENTE_CONVECCION_EXTERIOR,
                        T_e_max=self.tab_datos_entrada.parametros.TEMPERATURA_EXTERNA_MAXIMA, rango_temp=self.tab_datos_entrada.parametros.RANGO_TEMPERATURA,
                        T_i=self.tab_datos_entrada.parametros.TEMPERATURA_INTERIOR, zona_horaria=self.tab_datos_entrada.parametros.ZONA_HORARIA,
                        fecha=self.tab_datos_entrada.parametros.FECHA, latitud=self.tab_datos_entrada.parametros.LATITUD, 
                        orientacion=self.tab_datos_entrada.parametros.ORIENTACION, inclinacion=self.tab_datos_entrada.parametros.INCLINACION)

            # PERFIL DE TEMPERATURAS
            for line in self.tab_perfiles_temperatura.ax_temp_muro.lines[:]:
                line.remove()
            for line in self.tab_perfiles_temperatura.ax_temp_superficie.lines[:]:
                line.remove()
            
            # PARÁMETROS
            slider_temp_valor = int(self.tab_perfiles_temperatura.slider_simulacion.get())
            self.tab_perfiles_temperatura.min_temp = int(min(min(i) for i in self.resultados_simulacion[0]))
            self.tab_perfiles_temperatura.max_temp = int(max(max(i) for i in self.resultados_simulacion[0])) + 1
           
            # GRÁFICO DE LAS SUPERFICIES
            self.tab_perfiles_temperatura.ax_temp_superficie.plot(self.resultados_simulacion[-1], self.resultados_simulacion[1], color="C0", label="Superficie exterior")
            self.tab_perfiles_temperatura.ax_temp_superficie.plot(self.resultados_simulacion[-1], self.resultados_simulacion[2], color="C2", label="Superficie intermedia")
            self.tab_perfiles_temperatura.ax_temp_superficie.plot(self.resultados_simulacion[-1], self.resultados_simulacion[3], color="C3", label="Superficie interior")
            self.tab_perfiles_temperatura.ax_temp_superficie.legend()

            # GRÁFICO PERFIL DEL MURO
            self.tab_perfiles_temperatura.ax_temp_muro.plot(self.resultados_simulacion[-2], self.resultados_simulacion[0][slider_temp_valor], color="C1")
            self.tab_perfiles_temperatura.ax_temp_muro.set_ylim(self.tab_perfiles_temperatura.min_temp, self.tab_perfiles_temperatura.max_temp)
            self.tab_perfiles_temperatura.ax_temp_muro.set_yticks(np.linspace(self.tab_perfiles_temperatura.min_temp, self.tab_perfiles_temperatura.max_temp, 5))
            self.tab_perfiles_temperatura.ax_temp_muro.set_xlim(self.resultados_simulacion[-2][0],self.resultados_simulacion[-2][-1])
            acum_espesor = 0
            for i in self.tab_datos_entrada.parametros.ESPESORES:
                acum_espesor += i
                self.tab_perfiles_temperatura.ax_temp_muro.plot([(acum_espesor)*1000,(acum_espesor)*1000], [self.tab_perfiles_temperatura.min_temp, self.tab_perfiles_temperatura.max_temp], linestyle="dashed", lw=1, color="C3")
            
            
            self.tab_perfiles_temperatura.slider_simulacion.configure(state=customtkinter.NORMAL)
            self.tab_perfiles_temperatura.canvas.draw_idle()

            # FLUJOS DE CALOR
            for line in self.tab_resultados_termicos.ax_flujo_calor.lines[:]:
                line.remove()
            for line in self.tab_resultados_termicos.ax_ganancia_carga_termica.lines[:]:
                line.remove()

            # PARÁMETROS
            slider_flujo_valor = int(self.tab_resultados_termicos.slider_simulacion.get())
            self.tab_resultados_termicos.min_flujo = int(min(min(i) for i in self.resultados_simulacion[-3]))
            self.tab_resultados_termicos.max_flujo = int(max(max(i) for i in self.resultados_simulacion[-3])) + 1

            # GRÁFICO DE CARGA Y GANANCIA
            self.tab_resultados_termicos.ax_ganancia_carga_termica.plot(self.resultados_simulacion[-1], self.resultados_simulacion[-5], color="C0", label="Ganancia térmica")
            self.tab_resultados_termicos.ax_ganancia_carga_termica.plot(self.resultados_simulacion[-1], self.resultados_simulacion[-4], color="C2", label="Carga térmica")
            self.tab_resultados_termicos.ax_ganancia_carga_termica.legend()
            
            # GRÁFICO FLUJO CALOR
            self.tab_resultados_termicos.ax_flujo_calor.plot(self.resultados_simulacion[-2], self.resultados_simulacion[-3][slider_flujo_valor], color="C1")
            self.tab_resultados_termicos.ax_flujo_calor.set_ylim(self.tab_resultados_termicos.min_flujo, self.tab_resultados_termicos.max_flujo)
            self.tab_resultados_termicos.ax_flujo_calor.set_yticks(np.linspace(self.tab_resultados_termicos.min_flujo, self.tab_resultados_termicos.max_flujo, 5))
            self.tab_resultados_termicos.ax_flujo_calor.set_xlim(self.resultados_simulacion[-2][0], self.resultados_simulacion[-2][-1])
            
            acum_espesor = 0
            for i in self.tab_datos_entrada.parametros.ESPESORES:
                acum_espesor += i
                self.tab_resultados_termicos.ax_flujo_calor.plot([(acum_espesor)*1000,(acum_espesor)*1000],[self.tab_resultados_termicos.min_flujo,  self.tab_resultados_termicos.max_flujo], linestyle="dashed", lw=1, color="C3")
            self.tab_resultados_termicos.slider_simulacion.configure(state=customtkinter.NORMAL)
            self.tab_resultados_termicos.canvas.draw_idle()
            Alerta(self.master.master, text="Simulación númerica exitosa", time=1000, type_of_alert=0)
        except:
            Alerta(self.master.master, text="Algo salio mal, revisar los parámetros de entrada y reintentar la simulación", time=1000, type_of_alert=-1)

    def solucion_numerica_multiple(self):
        try:
            self.resultados_simulacion_multiple = []
            # OBTENER RESULTADOS DE SIMULACIÓN MULTIPLE
            for muro in self.tab_comparacion_multiple.muros:
                self.resultados_simulacion_multiple.append(solver(espesores=muro[2].ESPESORES, rho_cp=muro[2].PRODUCTO_RHO_CP, 
                            conductividades=muro[2].CONDUCTIVIDADES, difusividades=muro[2].DIFUSIVIDADES,
                            h_i=muro[2].COEFICIENTE_CONVECCION_INTERIOR, h_e=muro[2].COEFICIENTE_CONVECCION_EXTERIOR,
                            T_e_max=muro[2].TEMPERATURA_EXTERNA_MAXIMA, rango_temp=muro[2].RANGO_TEMPERATURA,
                            T_i=muro[2].TEMPERATURA_INTERIOR, zona_horaria=muro[2].ZONA_HORARIA,
                            fecha=muro[2].FECHA, latitud=muro[2].LATITUD, 
                            orientacion=muro[2].ORIENTACION, inclinacion=muro[2].INCLINACION))
            # GANANICA Y CARGA TÉRMICA
            for line in self.tab_comparacion_multiple.ax_ganancia_termica.lines[:]:
                line.remove()
            for line in self.tab_comparacion_multiple.ax_carga_termica.lines[:]:
                line.remove()

            for i in range(len(self.tab_comparacion_multiple.muros)):
                # GRÁFICO DE GANANCIA
                self.tab_comparacion_multiple.ax_ganancia_termica.plot(self.resultados_simulacion_multiple[i][-1], self.resultados_simulacion_multiple[i][-5], color=f"C{i}", label=self.tab_comparacion_multiple.muros[i][1])
                self.tab_comparacion_multiple.ax_ganancia_termica.legend()
                
                # GRÁFICO DE CARGA
                self.tab_comparacion_multiple.ax_carga_termica.plot(self.resultados_simulacion_multiple[i][-1], self.resultados_simulacion_multiple[i][-4], color=f"C{i}", linestyle="--", label=self.tab_comparacion_multiple.muros[i][1])
                self.tab_comparacion_multiple.ax_carga_termica.legend()

            self.tab_comparacion_multiple.canvas.draw_idle()
            Alerta(self.master.master, text="Simulación númerica exitosa", time=1000, type_of_alert=0)
        except:
            Alerta(self.master.master, text="Algo salio mal, revisar los parámetros de entrada y reintentar la simulación", time=1000, type_of_alert=-1)
# APLICACIÓN
class App(customtkinter.CTk):
    def __init__(self):
        super().__init__()
        
        self.grid_rowconfigure(0, weight=1)  # configure grid system
        self.grid_columnconfigure(0, weight=1)

        self.tab_view = MyTabView(master=self)
        self.tab_view.configure(corner_radius=15, segmented_button_selected_color="#083a5e", segmented_button_selected_hover_color="#022238")
        self.tab_view.grid(row=0, column=0, padx=20, pady=10, sticky="nsew")

customtkinter.set_appearance_mode("dark")
app = App()
>>>>>>> eedad3d08c27edbca3f97c278c7ca27170132576
app.mainloop()
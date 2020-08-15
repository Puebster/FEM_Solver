import numpy as np
import scipy
import math


class Knoten(object):
	"""docstring for Knoten"""
	def __init__(self,nummer, x, y):
		self.nummer = int(nummer)
		self.x = x
		self.y = y

class Element(object):
	"""docstring for Element"""
	def __init__(self, nummer, A_Knoten, E_Knoten, E, A):
		self.nummer = nummer
		self.anfang = A_Knoten
		self.ende = E_Knoten
		self.lange = math.sqrt(math.pow(self.ende.x - self.anfang.x, 2) + math.pow(self.ende.y - self.anfang.y, 2))
		self.alpha = math.atan2((self.ende.y - self.anfang.y),(self.ende.x - self.anfang.x))
		self.Steif_Matrix = (E*A/(self.lange))*np.matrix([[math.pow(math.cos(self.alpha),2),math.sin(self.alpha)*math.cos(self.alpha), -math.pow(math.cos(self.alpha),2), -math.sin(self.alpha)*math.cos(self.alpha)],[math.sin(self.alpha)*math.cos(self.alpha), math.pow(math.sin(self.alpha),2), -math.sin(self.alpha)*math.cos(self.alpha), -math.pow(math.sin(self.alpha),2)],[-math.pow(math.cos(self.alpha),2),-math.sin(self.alpha)*math.cos(self.alpha), math.pow(math.cos(self.alpha),2), math.sin(self.alpha)*math.cos(self.alpha)],[-math.sin(self.alpha)*math.cos(self.alpha), -math.pow(math.sin(self.alpha),2), math.sin(self.alpha)*math.cos(self.alpha), math.pow(math.sin(self.alpha),2)]])
		self.Elemt_trans_Matrix = np.matrix([[math.cos(self.alpha),math.sin(self.alpha),0,0],[0,0,math.cos(self.alpha),math.sin(self.alpha)]])
		self.C_matrix = np.matrix([[1, -1],[-1,1]])

	def calc_steif_mat(self):
		return self.Steif_Matrix

	def calc_Stabkrafte(self,verschiebunggsvektor):
		Stabkrafte = np.matmul(self.Steif_Matrix, verschiebunggsvektor)
		return Stabkrafte



def initiate_Knoten():
	Knoten_d = {}
	with open("Knotendatei.txt", "r") as Knotenliste:
		for	line in Knotenliste:
			daten = line.split()
			Knot = Knoten(daten[0],np.float64(daten[1]),np.float64(daten[2]))
			Knoten_d[daten[0]] = Knot
	Knotenliste.close()
	return Knoten_d

def initiate_Elemente(Knoten_d, E ,A):
	Elemente_d = {}
	with open("Elementdatei.txt", "r") as Elementliste:
		for line in Elementliste:
			daten = line.split()
			Eee = Element(daten[0], Knoten_d[daten[1]], Knoten_d[daten[2]], E , A)
			Elemente_d[daten[0]] = Eee
	Elementliste.close()
	return Elemente_d

def make_Zuordnungmatrix(Element_d, knotenanzahl):
	Zuordnung = []

	for anzahl in Element_d:
		Zuordnung_anfang_x = [0]*knotenanzahl*2
		Zuordnung_anfang_y = [0]*knotenanzahl*2
		Zuordnung_anfang_x[2*(Element_d[anzahl].anfang.nummer)-2] = 1
		Zuordnung_anfang_y[2*(Element_d[anzahl].anfang.nummer)-1] = 1
		Zuordnung.append(Zuordnung_anfang_x)
		Zuordnung.append(Zuordnung_anfang_y)

		Zuordnung_ende_x = [0]*knotenanzahl*2
		Zuordnung_ende_y = [0]*knotenanzahl*2
		Zuordnung_ende_x[2*(Element_d[anzahl].ende.nummer)-2] = 1
		Zuordnung_ende_y[2*(Element_d[anzahl].ende.nummer)-1] = 1
		Zuordnung.append(Zuordnung_ende_x)
		Zuordnung.append(Zuordnung_ende_y)

	Zuordnungsmatrix = np.matrix(Zuordnung)

	return Zuordnungsmatrix


def make_gesamt_steif_mat(Elemente_d, lenn):
	

	a = np.zeros((lenn*2,lenn*2))

	# print(Elemente_d["1"].Steif_Matrix)

	for i in Elemente_d:
		pos1 = int(Elemente_d[i].anfang.nummer) - 1
		pos2 = int(Elemente_d[i].ende.nummer) - 1

		a[pos1*2][pos1*2] += Elemente_d[i].Steif_Matrix.item((0,0))
		a[pos1*2][pos1*2+1] += Elemente_d[i].Steif_Matrix.item((0,1))
		a[pos1*2+1][pos1*2] += Elemente_d[i].Steif_Matrix.item((1,0))
		a[pos1*2+1][pos1*2+1] += Elemente_d[i].Steif_Matrix.item((1,1))

		a[pos2*2][pos1*2] += Elemente_d[i].Steif_Matrix.item((2,0))
		a[pos2*2][pos1*2+1] += Elemente_d[i].Steif_Matrix.item((2,1))
		a[pos2*2+1][pos1*2] += Elemente_d[i].Steif_Matrix.item((3,0))
		a[pos2*2+1][pos1*2+1] += Elemente_d[i].Steif_Matrix.item((3,1))


		a[pos1*2][pos2*2] += Elemente_d[i].Steif_Matrix.item((0,2))
		a[pos1*2][pos2*2+1] += Elemente_d[i].Steif_Matrix.item((0,3))
		a[pos1*2+1][pos2*2] += Elemente_d[i].Steif_Matrix.item((1,2))
		a[pos1*2+1][pos2*2+1] += Elemente_d[i].Steif_Matrix.item((1,3))


		a[pos2*2][pos2*2] += Elemente_d[i].Steif_Matrix.item((2,2))
		a[pos2*2][pos2*2+1] += Elemente_d[i].Steif_Matrix.item((2,3))
		a[pos2*2+1][pos2*2] += Elemente_d[i].Steif_Matrix.item((3,2))
		a[pos2*2+1][pos2*2+1] += Elemente_d[i].Steif_Matrix.item((3,3))

	return np.matrix(a)



A = 0.000314159265358979 # m Querschnittsfläche bei durchschnitt 2cm
E = 2000000000.0 

#Knotenlasten:
Knotenlasten = {}
Knotenlasten["P1_x"] = 50000. 
Knotenlasten["P12_x"] = -40000. 
Knotenlasten["P12_y"] = -60000. 
Knotenlasten["P17_y"] = 40000. 


#Geometrische Randbedingunen
geom_rand = {}
geom_rand["v10_x"] = 0.
geom_rand["v10_y"] = 0.
geom_rand["v18_x"] = 0.
geom_rand["v18_y"] = 0.

Knots = initiate_Knoten()
Elem = initiate_Elemente(Knots, E, A)


######
v = []
p = []
for knotenv in range(len(Knots)*2):
	v.append([-1])
	p.append([0])
indices_rand = []
# v = [-1] * len(Knots) * 2
for randbeding in geom_rand:
	strin = str(randbeding)
	check = strin.split("_")
	check[0] = int(check[0].replace("v", ""))
	if check[1] == "x":
		check[0] = (check[0] - 1) * 2
	else:
		check[0] = (check[0] - 1) * 2 + 1
	indices_rand.append(check[0])
	v[check[0]][0] = geom_rand[randbeding]
	p[check[0]][0] = "a"
indices_rand_c = list(range(len(Knots)*2))

for kn_last in Knotenlasten:
	strin = str(kn_last)
	check = strin.split("_")
	check[0] = int(check[0].replace("P", ""))
	if check[1] == "x":
		check[0] = (check[0] - 1) * 2
	else:
		check[0] = (check[0] - 1) * 2 + 1
	p[check[0]][0] = Knotenlasten[kn_last]

for j in indices_rand:
	indices_rand_c.remove(j)
######

Zuordnungsmatrix = make_Zuordnungmatrix(Elem,len(Knots))
f_steif_mat = make_gesamt_steif_mat(Elem, len(Knots))

K1 = f_steif_mat.take(indices_rand_c, axis = 0).take(indices_rand_c,axis = 1)

p_zerlegt =  np.array(p).take(indices_rand_c)
p_zerlegt = np.transpose([p_zerlegt]).astype(np.float64)

rest_of_v = np.linalg.solve(np.asarray(K1), p_zerlegt)
np.set_printoptions(suppress=True)
#rest_of_v = groß_K_zerlegt_1.getI() * np.matrix(p)[indices_rand_c]

K2 = f_steif_mat.take(indices_rand, axis = 0).take(indices_rand_c,axis = 1)
groß_K_zerlegt_2 = f_steif_mat[indices_rand]
groß_K_zerlegt_2 = f_steif_mat[:,indices_rand_c]
rest_of_p = K2 * rest_of_v


einfacher_counter = 0
for eintrag in range(len(v)):
	if v[eintrag][0] == -1:
		v[eintrag][0] = float(rest_of_v[einfacher_counter][0])
		einfacher_counter += 1

einfacher_counter = 0
for eintrag in range(len(p)):
	if p[eintrag][0] == "a":
		p[eintrag][0] = float(rest_of_p[einfacher_counter][0])
		einfacher_counter += 1

with open("results_model_2v1.txt", "a") as output:

	output.write("Knotenverschiebungsvektor v:\n\n")
	output.write(str(np.matrix(v)))
	output.write("\n\n")

	print("Knotenverschiebungsvektor v:")
	print(np.matrix(v))
	print("\n")
	print("\n")


	output.write("Lagerreaktionskräfte p:\n\n")
	output.write(str(np.matrix(p)))
	output.write("\n\n")
	print("Lagerreaktionskräfte p:")
	print(np.matrix(p))
	print("\n")
	print("\n")

	output.write("Stabdehnung:" + "\t\t" + "Stabspannung\n\n")
	for letzte in Elem:
		Stabk = (((E*A)/Elem[letzte].lange)*Elem[letzte].C_matrix) * Elem[letzte].Elemt_trans_Matrix * Zuordnungsmatrix[[(int(letzte)-1)*4,(int(letzte)-1)*4+1,(int(letzte)-1)*4+2,(int(letzte)-1)*4+3]]*np.matrix(v)
		# output.write("\nStabspannung " + letzte + ":\n\t")
		# output.write(str(Stabk[1]/A).replace("]", "").replace("[", ""))
		# output.write("\nStabdehnung " + letzte + ":\n\t")
		# output.write(str((Stabk[1]/A)/E).replace("]", "").replace("[", ""))

		output.write(str((Stabk[1]/A)/E).replace("]", "").replace("[", "") + "\t" + "\t" + str(Stabk[1]/A).replace("]", "").replace("[", "")+ "\n")


		print("Stabkraft von Stab " + letzte + ":")
		print(Stabk[1]/Elem[letzte].lange)
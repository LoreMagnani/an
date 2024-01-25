# usare nella stessa cartella in cui ci sono i file lhe. Questi devono chiamarsi nome_numero.lhe. Prende in argomento il nome dei file lhe.
# Legge i vari file lhe, calcola le variabili interessate, salva i valori in un file .root che è un TTree 
import sys
sys.path.insert(0, '/..')
import ROOT
import pylhe
import awkward as ak
import numpy as np
import os
import vector
vector.register_awkward()

#opera una selezione degli eventi
def selection(arr , i):
    if  arr[i].particles[-3].vector.pt > 20 \
            and arr[i].particles[-4].vector.pt > 20 \
            and abs(arr[i].particles[-4].vector.eta) <= 2.4 \
            and abs(arr[i].particles[-3].vector.eta) <= 2.4 \
            and arr[i].particles[-2].vector.pt >= 30 \
            and arr[i].particles[-1].vector.pt >= 30 \
            and abs(arr[i].particles[-1].vector.eta) <= 4.7 \
            and abs(arr[i].particles[-2].vector.eta) <= 4.7 \
            and (arr[i].particles[-4].vector+arr[i].particles[-3].vector).M > 50 \
            and (arr[i].particles[-2].vector+arr[i].particles[-1].vector).M > 200:
        if arr[i].particles[-3].id == 11 and arr[i].particles[-4].id == -11:
            return True
        elif arr[i].particles[-3].id == 13 and arr[i].particles[-4].id == -13:
            return True
    else:
        return False

def main():
    #ottiene i percorsi assoluti dei file lhe e li salva in un array
    folder_path = os.getcwd()
    prefix = sys.argv[1]
    extension = ".lhe"
    unsorted_files = [file for file in os.listdir(folder_path) if file.startswith(prefix) and file.endswith(extension)]
    unsorted_files = np.array(unsorted_files)
    def custom_sort(file_name):
        return int(file_name.split('_')[1].split('.')[0])
    files = sorted(unsorted_files, key=custom_sort)
    files = list(map(lambda k: folder_path + "/" + k, files))

    #questo array contiene le sezioni d'urto di ogni file
    xs = []

    #itero su ogni file contenuto in files
    for index in range(len(files)):
        file = files[index]
 
        print(file, index)

        #salvo la sezione d'urto
        xs.append(pylhe.read_lhe_init(file)['procInfo'][0]['xSection'])

        #trova nomi dei pesi
        weight_names = ['rwgt_' + str(i) for i in range(1, 46)]

        # Lettura dei dati LHE
        arr = pylhe.to_awkward(pylhe.read_lhe_with_attributes(file) , weight_names)

        #calcolo peso nominale per tutti gli eventi e controllo la somma restituisca la sezione d'urto
        tot_nom_weights = arr.eventinfo.weight * xs[index] * 1000/ ak.sum(arr.eventinfo.weight)
        print("cross-section del file: " + str(xs[index]))
        print("somma pesi nominali: " + str(ak.sum(tot_nom_weights)))

        # Creazione di un file ROOT e di un TTree
        output_file = ROOT.TFile("output_" + str(index) + ".root" , "RECREATE")
        output_tree = ROOT.TTree("events", "Event tree")

        # Dichiarazione delle variabili da inserire nel TTree
        mll = np.zeros(1, dtype=float)
        deltaPhill = np.zeros(1, dtype=float)
        deltaEtall = np.zeros(1, dtype=float)
        ptl1 = np.zeros(1, dtype=float)
        phil1 = np.zeros(1, dtype=float)
        etal1 = np.zeros(1, dtype=float)
        ptl2 = np.zeros(1, dtype=float)
        phil2 = np.zeros(1, dtype=float)
        etal2 = np.zeros(1, dtype=float)
        mJJ = np.zeros(1, dtype=float)
        deltaPhiJJ = np.zeros(1, dtype=float)
        deltaEtaJJ = np.zeros(1, dtype=float)
        ptJ1 = np.zeros(1, dtype=float)
        phiJ1 = np.zeros(1, dtype=float)
        etaJ1 = np.zeros(1, dtype=float)
        ptJ2 = np.zeros(1, dtype=float)
        phiJ2 = np.zeros(1, dtype=float)
        etaJ2 = np.zeros(1, dtype=float)
        weights = np.zeros(45 , dtype=float)
        nom_weights = np.zeros(1 , dtype=float)

        # Creazione delle branch nel TTree
        output_tree.Branch("mll" , mll , "mll/D")
        output_tree.Branch("deltaPhill" , deltaPhill , "deltaPhill/D")
        output_tree.Branch("deltaEtall" , deltaEtall , "deltaEtall/D")
        output_tree.Branch("ptl1" , ptl1 , "ptl1/D")
        output_tree.Branch("phil1" , phil1 , "phil1/D")
        output_tree.Branch("etal1" , etal1 , "etal1/D")
        output_tree.Branch("ptl2" , ptl2 , "ptl2/D")
        output_tree.Branch("phil2" , phil2 , "phil2/D")
        output_tree.Branch("etal2" , etal2 , "etal2/D")
        output_tree.Branch("mJJ" , mJJ , "mJJ/D")
        output_tree.Branch("deltaPhiJJ" , deltaPhiJJ , "deltaPhiJJ/D")
        output_tree.Branch("deltaEtaJJ" , deltaEtaJJ , "deltaEtaJJ/D") 
        output_tree.Branch("ptJ1" , ptJ1 , "ptJ1/D")
        output_tree.Branch("phiJ1" , phiJ1 , "phiJ1/D")
        output_tree.Branch("etaJ1" , etaJ1 , "etaJ1/D") 
        output_tree.Branch("ptJ2" , ptJ2 , "ptJ2/D")
        output_tree.Branch("phiJ2" , phiJ2 , "phiJ2/D")
        output_tree.Branch("etaJ2" , etaJ2 , "etaJ2/D")
        output_tree.Branch("weights" , weights , "weights[45]/D")
        output_tree.Branch("nom_weights" , nom_weights , "nom_weights/D")
 
        #itero sugli eventi
        for i in range(len(arr)):
            if selection(arr , i):

                #-------------------Leptons variables-----------------------------

                mll[0] = (arr[i].particles[-4].vector+arr[i].particles[-3].vector).M
                
                deltaPhill[0] = (arr[i].particles[-4].vector).deltaphi(arr[i].particles[-3].vector)
                deltaEtall[0] = (arr[i].particles[-4].vector).deltaeta(arr[i].particles[-3].vector)

                firstIsFirst_l = np.array(arr[i].particles[-4].vector.pt > arr[i].particles[-3].vector.pt)

                ptl1[0]  = np.array(arr[i].particles[-4].vector.pt)
                phil1[0] = np.array(arr[i].particles[-4].vector.phi)
                etal1[0] = np.array(arr[i].particles[-4].vector.eta)
                ptl1[0] = np.where(firstIsFirst_l, arr[i].particles[-4].vector.pt, arr[i].particles[-3].vector.pt)
                phil1[0] = np.where(firstIsFirst_l, arr[i].particles[-4].vector.phi, arr[i].particles[-3].vector.phi)
                etal1[0] = np.where(firstIsFirst_l, arr[i].particles[-4].vector.eta, arr[i].particles[-3].vector.eta)

                
                #ptl1[0] [~firstIsFirst_l] = np.array(arr[i].particles[-3].vector.pt) [~firstIsFirst_l]
                #phil1[0][~firstIsFirst_l] = np.array(arr[i].particles[-3].vector.phi)[~firstIsFirst_l]
                #etal1[0][~firstIsFirst_l] = np.array(arr[i].particles[-3].vector.eta)[~firstIsFirst_l]

                ptl2[0]  = np.array(arr[i].particles[-3].vector.pt)
                phil2[0] = np.array(arr[i].particles[-3].vector.phi)
                etal2[0] = np.array(arr[i].particles[-3].vector.eta)
                ptl2[0] = np.where(firstIsFirst_l, arr[i].particles[-3].vector.pt, arr[i].particles[-4].vector.pt)
                phil2[0] = np.where(firstIsFirst_l, arr[i].particles[-3].vector.phi, arr[i].particles[-4].vector.phi)
                etal2[0] = np.where(firstIsFirst_l, arr[i].particles[-3].vector.eta, arr[i].particles[-4].vector.eta)


                #ptl2[0] [~firstIsFirst_l] = np.array(arr[i].particles[-4].vector.pt) [~firstIsFirst_l]
                #phil2[0][~firstIsFirst_l] = np.array(arr[i].particles[-4].vector.phi)[~firstIsFirst_l]
                #etal2[0][~firstIsFirst_l] = np.array(arr[i].particles[-4].vector.eta)[~firstIsFirst_l]

                #-------------------Jet variables-------------------------------

                mJJ[0] = (arr[i].particles[-2].vector+arr[i].particles[-1].vector).M

                firstIsFirst_J = np.array(arr[i].particles[-2].vector.pt > arr[i].particles[-1].vector.pt)

                deltaPhiJJ[0] = (arr[i].particles[-2].vector).deltaphi(arr[i].particles[-1].vector)
                deltaEtaJJ[0] = (arr[i].particles[-2].vector).deltaeta(arr[i].particles[-1].vector)

                ptJ1[0] = np.array(arr[i].particles[-2].vector.pt)
                phiJ1[0] = np.array(arr[i].particles[-2].vector.phi)
                etaJ1[0] = np.array(arr[i].particles[-2].vector.eta)
                ptJ1[0] = np.where(firstIsFirst_J, arr[i].particles[-2].vector.pt, arr[i].particles[-1].vector.pt)
                phiJ1[0] = np.where(firstIsFirst_J, arr[i].particles[-2].vector.phi, arr[i].particles[-1].vector.phi)
                etaJ1[0] = np.where(firstIsFirst_J, arr[i].particles[-2].vector.eta, arr[i].particles[-1].vector.eta)

                #ptJ1[0] [~firstIsFirst_J] =np.array( arr[i].particles[-1].vector.pt) [~firstIsFirst_J]
                #phiJ1[0] [~firstIsFirst_J] = np.array(arr[i].particles[-1].vector.phi) [~firstIsFirst_J]
                #etaJ1[0] [~firstIsFirst_J] = np.array(arr[i].particles[-1].vector.eta) [~firstIsFirst_J]

                ptJ2[0] =np.array(arr[i].particles[-1].vector.pt)
                phiJ2[0]  =np.array(arr[i].particles[-1].vector.phi)
                etaJ2[0] = np.array(arr[i].particles[-1].vector.eta)
                ptl1[0] = np.where(firstIsFirst_l, arr[i].particles[-1].vector.pt, arr[i].particles[-2].vector.pt)
                phiJ2[0] = np.where(firstIsFirst_l, arr[i].particles[-1].vector.phi, arr[i].particles[-2].vector.phi)
                etaJ2[0] = np.where(firstIsFirst_l, arr[i].particles[-1].vector.eta, arr[i].particles[-2].vector.eta)

                #ptJ2[0] [~firstIsFirst_J] =np.array( arr[i].particles[-2].vector.pt) [~firstIsFirst_J]
                #phiJ2[0] [~firstIsFirst_J] = np.array(arr[i].particles[-2].vector.phi) [~firstIsFirst_J]
                #etaJ2[0] [~firstIsFirst_J] = np.array(arr[i].particles[-2].vector.eta) [~firstIsFirst_J]

                #-------------------weights variables-------------------------------

                weights[:] = arr[i].weights.values
                nom_weights[0] = tot_nom_weights[i]

                output_tree.Fill()

        output_file.Write()
        output_file.Close()

    # Creazione di un file ROOT e di un TTree per la sezione d'urto
    output_file_xs = ROOT.TFile("sezione_urto.root", "RECREATE")
    output_tree_xs = ROOT.TTree("xs_tree", "Sezione d'urto tree")

    # Dichiarazione della variabile da inserire nel TTree
    xs_value = np.zeros(1, dtype=float)

    # Creazione della branch nel TTree
    output_tree_xs.Branch("xs_value", xs_value, "xs_value/D")

    # Riempimento del TTree con i valori della sezione d'urto
    for value in xs:
        xs_value[0] = value
        output_tree_xs.Fill()

    # Scrittura su file e chiusura del file ROOT per la sezione d'urto
    output_file_xs.Write()
    output_file_xs.Close()

    print(xs)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Funzionamento: Python3 an.py -nome_dei_file-")
        sys.exit(1)

    main()






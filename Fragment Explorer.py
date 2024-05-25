import sys
import os
import csv
import scaffoldgraph
from rdkit import Chem
from rdkit.Chem import SDMolSupplier
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from rdkit.Chem import Draw, Descriptors, Lipinski, Recap, BRICS, rdFMCS, AllChem
import matplotlib.pyplot as plt
from PyQt5.QtWidgets import QMessageBox, QApplication, QMainWindow, QStatusBar, QPushButton, QAbstractItemView, QWidget, QVBoxLayout, QFileDialog, QAction, QTextEdit, QCheckBox, QSplitter, QGroupBox, QLabel, QTableWidget, QTableWidgetItem, QScrollArea, QRadioButton,QInputDialog
from PyQt5.QtGui import QImage, QPixmap, QFont
from PyQt5.QtCore import Qt, QUrl
from PyQt5.QtNetwork import QNetworkAccessManager, QNetworkRequest
from io import BytesIO
import pandas as pd
from scaffoldgraph import ScaffoldNetwork
from PyQt5.QtCore import QTimer
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdFMCS
import os
from PyQt5.QtGui import QPixmap
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QPushButton, QLabel, QMessageBox
from io import BytesIO
from PyQt5.QtCore import Qt

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdFMCS
from rdkit.Chem import AllChem
import os
from PyQt5.QtGui import QPixmap
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QPushButton, QLabel, QMessageBox
from io import BytesIO
from PyQt5.QtCore import Qt

from rdkit.Chem import SDMolSupplier, MolToSmiles, MolFromSmiles
from rdkit.Chem import Draw
from PyQt5.QtGui import QPixmap, QImage
from PyQt5.QtCore import Qt

class Tfrag(QMainWindow):
   
    def __init__(self):
        super().__init__()
        self.splitter = None
        self.left_frame = None
        #self.center_frame = None
        self.right_frame = None
        self.initUI()
        self.features = QTableWidget()
        self.selected_fragments = set()
        self.brics_fragments = set()
        self.recap_fragments = set()
        self.selected_fragment = ""

   
    
    
    def displaySDF(self, sdf_file):
     # Read the SDF file and visualize the first molecule
     suppl = Chem.SDMolSupplier(sdf_file)
     for mol in suppl:
        if mol is not None:
            # Generate a 2D image of the molecule with a fixed size of 300x300 pixels
            img = Draw.MolToImage(mol, size=(300, 300))
            buffer = BytesIO()
            img.save(buffer, format="PNG")
            # Convert image to QPixmap and display it
            qt_pixmap = QPixmap()
            qt_pixmap.loadFromData(buffer.getvalue())
            # Set the pixmap to the label without scaling
            self.image_display.setPixmap(qt_pixmap)
            self.image_display.setFixedSize(300, 300)  # Set a fixed size for the label
            break  # Only display the first molecule for now
       
   
    def openFileDialog(self):
     # Open file dialog and get the selected file path, filtering for SDF files only
     options = QFileDialog.Options()
     fileName, _ = QFileDialog.getOpenFileName(self, "Select SDF File", "", "SDF Files (*.sdf)", options=options)

     if fileName:
        print("Selected file:", fileName)
        # Convert SDF file to SMILES strings
        smiles_list = self.convertSDFtoSMILES(fileName)
        # Display the SMILES or fragments in the table widget
        self.displaySMILES(smiles_list)
        # Update the molecule visualization
        self.updateMoleculeFromSDF(smiles_list)
     
    
    def openFileDialog(self):
     # Open file dialog and get the selected file path, filtering for SDF files only
     options = QFileDialog.Options()
     fileName, _ = QFileDialog.getOpenFileName(self, "Select SDF File", "", "SDF Files (*.sdf)", options=options)

     if fileName:
        print("Selected file:", fileName)
        # Convert SDF file to SMILES strings
        smiles_list = self.convertSDFtoSMILES(fileName)
        # Display the SMILES or fragments in the table widget
        self.displaySMILES(smiles_list)
        # Update the molecule visualization
        self.updateMoleculeFromSDF(smiles_list)
     
    
    def convertSDFtoSMILES(self, fileName):
        smiles_list = []
        suppl = SDMolSupplier(fileName)
        for mol in suppl:
            if mol is not None:
                smiles = MolToSmiles(mol)
                smiles_list.append(smiles)

                # Display image if smiles is available
                molecule_image = self.loadMolecule(smiles)
                if molecule_image:
                    self.displayMoleculeImage(molecule_image)

        return smiles_list
        
    def displayMoleculeImage(self, molecule_image):
        if molecule_image:
            self.image_display.setPixmap(molecule_image.scaled(self.image_display.size(), Qt.KeepAspectRatio))
        
    
   
    def displaySMILES(self, smiles_list):
     self.smiles_list = smiles_list  # Store the SMILES list for later use
     self.table_widget.setRowCount(len(smiles_list))
     self.table_widget.setColumnCount(1)
     self.table_widget.setHorizontalHeaderLabels(["SMILES"])

     for row, smiles in enumerate(smiles_list):
        item = QTableWidgetItem(smiles)
        self.table_widget.setItem(row, 0, item)  # Display each SMILES in the table widget

     # Increase the column width to ensure full view of SMILES
     self.table_widget.setColumnWidth(0, 1500)  # Adjust the width (e.g., 300 pixels) as needed


   
   
    def updateMoleculeFromSDF(self, smiles_list):
     if smiles_list:
        # For simplicity, take the first SMILES from the list
        smiles = smiles_list[0]
        molecule_image = self.loadMolecule(smiles)
        if molecule_image:
            self.image_display.setPixmap(molecule_image.scaled(self.image_display.size(), Qt.KeepAspectRatio))
   
    def displaySDF(self, sdf_file):
        # Read the SDF file and visualize the first molecule
        suppl = Chem.SDMolSupplier(sdf_file)
        for mol in suppl:
            if mol is not None:
                # Generate a 2D image of the molecule
                # Generate a 2D image of the molecule
                img = Draw.MolToImage(mol, size=(200, 200))  # Specify the image size here
               
                buffer = BytesIO()
                img.save(buffer, format="PNG")
                # Convert image to QPixmap and display it
                qt_pixmap = QPixmap()
                qt_pixmap.loadFromData(buffer.getvalue())
                self.image_display.setPixmap(qt_pixmap.scaled(self.image_display.size(), Qt.KeepAspectRatio))
                break  # Only display the first molecule for now
   


        
   
    def createLeftFrame(self):
        
      
        
        # Create left frame with text box and legends
        self.left_frame = QWidget()
        self.left_layout = QVBoxLayout(self.left_frame)
        self.left_frame.setFixedWidth(320)  # Set fixed width for the
        self.smiles_label = QLabel("SMILES String:")
        self.smiles_textbox = QTextEdit()
        self.smiles_textbox.setFixedHeight(140)
        self.left_layout.addWidget(self.smiles_label)
        self.left_layout.addWidget(self.smiles_textbox)
     
        self.smiles_textbox.textChanged.connect(self.updateMolecule)
       
        # Create a button to open the file dialog
        self.file_button = QPushButton("Choose File")
        self.left_layout.addWidget(self.file_button)
        self.file_button.clicked.connect(self.openFileDialog)


        # Fragment Splitter GroupBox
        self.fragment_splitter_box = QGroupBox("Fragment Splitter")
        self.fragment_splitter_layout = QVBoxLayout(self.fragment_splitter_box)
         
        # Button for RECAP fragmentation
        self.fragment_splitter_button_recap = QPushButton('Get fragments by RECAP')
        self.fragment_splitter_button_recap.clicked.connect(self.fragrecap)
        self.fragment_splitter_layout.addWidget(self.fragment_splitter_button_recap)
        
        # Button for BRICS fragmentation
        self.fragment_splitter_button_brics = QPushButton('Get fragments by BRICS')
        self.fragment_splitter_button_brics.clicked.connect(self.fragbrics)
        self.fragment_splitter_layout.addWidget(self.fragment_splitter_button_brics)
       
        self.fragment_splitter_button_compar = QPushButton('Download')
        self.fragment_splitter_button_compar.clicked.connect(self.compar)
        self.fragment_splitter_layout.addWidget(self.fragment_splitter_button_compar)
       
        self.left_layout.addWidget(self.fragment_splitter_box)
            
        # Legend Box 2: Fragment Filters
        self.fragment_filters_box = QGroupBox("Fragment Filters")
        
        self.fragment_filters_layout = QVBoxLayout(self.fragment_filters_box)
        self.filters = {
            "Rule of 3 Filter": self.filter_rule_of_3,
              }
        # Add buttons for each filter
        for filter_name in self.filters.keys():
             button = QPushButton(filter_name)
             button.clicked.connect(lambda checked, name=filter_name: self.apply_filter(name))
             self.fragment_filters_layout.addWidget(button)
        # Add this to your setupUi or __init__ method where you create UI components
        
        # Add the button to your layout or toolbar as needed

        self.left_layout.addWidget(self.fragment_filters_box)
        
       
                # Create the Similarity search box
        # Assuming this code is within a method of your QMainWindow or QWidget subclass

        self.similarity_search_box = QGroupBox("Scaffold Hopping")
        self.similarity_search_layout = QVBoxLayout(self.similarity_search_box)

        # Create the 'Scaffold Hopping' button and add it to the group box layout
        self.pubchem_button = QPushButton('Scaffold Hopping')
        self.pubchem_button.clicked.connect(self.scaffold_hopping)
        self.similarity_search_layout.addWidget(self.pubchem_button)

        # Create the 'Download Scaffold' button and add it to the group box layout
        self.download_button = QPushButton("Download Scaffold")
        self.download_button.clicked.connect(self.download_scaffolds)
        self.similarity_search_layout.addWidget(self.download_button)

        # Create the 'Download All Scaffolds' button and add it to the group box layout
        self.download_all_button = QPushButton("Download All Scaffolds")
        self.download_all_button.clicked.connect(self.download_all_scaffolds)
        self.similarity_search_layout.addWidget(self.download_all_button)

        # Add the group box to the left layout
        self.left_layout.addWidget(self.similarity_search_box)
        
        
        # Add the image display at the bottom of the left frame
        self.image_display = QLabel()
        self.image_display.setAlignment(Qt.AlignCenter)
        self.image_display.setFixedSize(300, 300)  # Set a fixed size for the image display
        self.left_layout.addWidget(self.image_display)
        
        # Set the central widget to the left frame which contains the left layout
        self.setCentralWidget(self.left_frame)  # Assuming self.left_frame is already defined
        font = QFont()
        font.setPointSize(10)  # Set the font size to 12
        for label in self.left_frame.findChildren(QLabel):
              label.setFont(font)        
        
    
    
    def download_all_scaffolds(self):
        all_scaffolds = []
        # Collect all scaffolds
        for row in range(self.table_widget.rowCount()):
            smiles = self.table_widget.item(row, 1).text()  # Assuming the SMILES are in the third column
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                all_scaffolds.append(mol)

        # Write all scaffolds to an SDF file
        with Chem.SDWriter('all_scaffolds.sdf') as writer:
            for mol in all_scaffolds:
                writer.write(mol)

        # Inform the user that the download is complete
        QMessageBox.information(self, "Download Complete", "All scaffolds have been downloaded as 'all_scaffolds.sdf'.")
    
    def download_scaffolds(self):
     selected_scaffolds = []
     # Collect the selected scaffolds
     for row in range(self.table_widget.rowCount()):
        checkbox = self.table_widget.cellWidget(row, 1)  # Assuming the checkbox is in the second column
        if checkbox is not None and checkbox.isChecked():
            smiles = self.table_widget.item(row, 2).text()  # Assuming the SMILES are in the third column
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                selected_scaffolds.append(mol)

     # Write the selected scaffolds to an SDF file
     with Chem.SDWriter('selected_scaffolds.sdf') as writer:
        for mol in selected_scaffolds:
            writer.write(mol)

     # Inform the user that the download is complete
     QMessageBox.information(self, "Download Complete", "Selected scaffolds have been downloaded as 'selected_scaffolds.sdf'.")


    def scaffold_hopping(self):
     selected_smiles = []

     # Collect SMILES strings of selected molecules
     for row in range(self.table_widget.rowCount()):
        checkbox_item = self.table_widget.cellWidget(row, 0)
        if checkbox_item and isinstance(checkbox_item, QCheckBox) and checkbox_item.isChecked():
            smiles_item = self.table_widget.item(row, 2)
            if smiles_item:
                smiles = smiles_item.text()
                selected_smiles.append(smiles)

     # Check if at least one fragment is selected
     if not selected_smiles:
        choice = QMessageBox.question(self, "No Fragments Selected",
                                       "No fragments selected. Do you want to input a molecule directly as SMILES - Yes or upload an SDF file - No?",
                                       QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel)
        if choice == QMessageBox.Yes:
            molecule, ok = QInputDialog.getText(self, "Input Molecule", "Enter the SMILES string:")
            if ok and molecule.strip():
                selected_smiles.append(molecule.strip())
                # Display image if SMILES is provided
                self.displayMoleculeImage(self.loadMolecule(molecule.strip()))
            else:
                QMessageBox.warning(self, "Invalid Input", "Please enter a valid SMILES string.")
                return
        elif choice == QMessageBox.No:
            options = QFileDialog.Options()
            fileName, _ = QFileDialog.getOpenFileName(self, "Select SDF File", "", "SDF Files (*.sdf)", options=options)
            if fileName:
                selected_smiles.extend(self.convertSDFtoSMILES(fileName))
            else:
                QMessageBox.warning(self, "Invalid Input", "Please select a valid SDF file.")
                return
        else:
            return

     # Generate scaffolds using ScaffoldGraph
     scaffolds = []
     for smiles in selected_smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            frags = scaffoldgraph.get_all_murcko_fragments(mol)
            scaffolds.extend(frags)

     # Display the filtered scaffolds in the right frame
     self.display_scaffolds(scaffolds)




    def display_scaffolds(self, scaffolds):
     # Clear existing data in the QTableWidget
     self.table_widget.setRowCount(0)

     # Define the column titles without the "No" column
     column_titles = ["Select", "SMILES", "Structure"]
     self.table_widget.setColumnCount(len(column_titles))
     self.table_widget.setHorizontalHeaderLabels(column_titles)

     # Define the desired height for the table rows to accommodate the new image size
     row_height = 300  # Set to the same height as the image

     # Populate QTableWidget with scaffold information
     for idx, scaffold_mol in enumerate(scaffolds):
        rowPosition = self.table_widget.rowCount()
        self.table_widget.insertRow(rowPosition)

        # Add checkbox for the scaffold in the first column
        checkbox = QCheckBox()
        checkbox.setChecked(False)
        self.table_widget.setCellWidget(rowPosition, 0, checkbox)

        # Convert the Mol object to a SMILES string
        scaffold_smiles = Chem.MolToSmiles(scaffold_mol)
        # Add scaffold SMILES
        self.table_widget.setItem(rowPosition, 1, QTableWidgetItem(scaffold_smiles))

        # Create a QLabel to display the scaffold image with a size of 300x300 pixels
        label = QLabel()
        if scaffold_mol:
            mol_image = Draw.MolToImage(scaffold_mol, size=(300, 300))
            qimage = QImage(mol_image.tobytes(), mol_image.width, mol_image.height, QImage.Format_RGB888)
            pixmap = QPixmap.fromImage(qimage)
            label.setPixmap(pixmap)
            self.table_widget.setCellWidget(rowPosition, 2, label)  # Set the widget in the cell

        # Set the height of the row to accommodate the new image size
        self.table_widget.setRowHeight(rowPosition, row_height)

     # Adjust column widths as needed
     self.table_widget.setColumnWidth(0, 50)  # Select column width
     self.table_widget.setColumnWidth(1, 200) # SMILES column width
     self.table_widget.setColumnWidth(2, 300) # Structure column width to match image width

     if not hasattr(self, 'similarity_mcs_button'):
        self.similarity_mcs_button = QPushButton('Similarity MCS')
        self.similarity_mcs_button.clicked.connect(self.calculate_mcs)
        self.right_layout.addWidget(self.similarity_mcs_button) 
    
    
    def SmilesMCStoGridImage(self, smiles, align_substructure=True, verbose=False, **kwargs):
        mols = [Chem.MolFromSmiles(smile) for smile in smiles]
        res = rdFMCS.FindMCS(mols, **kwargs)
        mcs_smarts = res.smartsString
        mcs_mol = Chem.MolFromSmarts(mcs_smarts)
        smarts_and_mols = [mcs_mol] + mols

        smarts_legend = "Max. substructure match"
        if isinstance(smiles, dict):
            mol_legends = [smiles[molecule] for molecule in smiles]
        else:
            mol_legends = ["" for mol in mols]

        legends = [smarts_legend] + mol_legends

        # Ensure each match is within the range of the molecule's atom indices
        matches = []
        for mol in mols:
            match = mol.GetSubstructMatch(mcs_mol)
            if all(idx < mol.GetNumAtoms() for idx in match):
                matches.append(match)
            else:
                matches.append([])  # Provide an empty match if indices are out of range

        subms = [x for x in smarts_and_mols if x.HasSubstructMatch(mcs_mol)]
        AllChem.Compute2DCoords(mcs_mol)

        if align_substructure:
            for m in subms:
                Draw.MolToImage(m, highlightAtoms=m.GetSubstructMatch(mcs_mol))

        # Ensure matches list has the correct length and valid atom indices
        matches = [[]] + matches

        drawing = Draw.MolsToGridImage(smarts_and_mols, highlightAtomLists=matches, legends=legends)

        if verbose:
            return drawing, mcs_smarts, mcs_mol, mols
        else:
            return drawing

    def calculate_mcs(self):
        selected_smiles = []
        for row in range(self.table_widget.rowCount()):
            checkbox = self.table_widget.cellWidget(row, 0)
            if checkbox.isChecked():
                smiles = self.table_widget.item(row, 1).text()
                selected_smiles.append(smiles)

        if len(selected_smiles) >= 2:
            grid_image, mcs_smarts, mcs_mol, mols = self.SmilesMCStoGridImage(selected_smiles, align_substructure=True, verbose=True)
            self.display_mcs_result(grid_image)
        else:
            QMessageBox.warning(self, "Selection Error", "Please select at least two molecules for MCS calculation.")

    def display_mcs_result(self, grid_image):
        if hasattr(self, 'mcs_result_container') and self.mcs_result_container is not None:
            self.right_layout.removeWidget(self.mcs_result_container)
            self.mcs_result_container.deleteLater()
            self.mcs_result_container = None

        self.mcs_result_container = QWidget()
        mcs_result_layout = QVBoxLayout(self.mcs_result_container)

        close_button = QPushButton('X')
        close_button.setFixedSize(20, 20)
        close_button.clicked.connect(self.close_mcs_result)
        mcs_result_layout.addWidget(close_button)

        image_bytes = BytesIO()
        grid_image.save(image_bytes, format='PNG')
        qt_pixmap = QPixmap()
        qt_pixmap.loadFromData(image_bytes.getvalue())
        
        image_path = os.path.join(os.path.expanduser('~'), 'Downloads', 'mcs_result.png')
        grid_image.save(image_path, 'PNG')

        image_label = QLabel()
        image_label.setPixmap(qt_pixmap.scaled(600, 1000, Qt.KeepAspectRatio))
        mcs_result_layout.addWidget(image_label)
        
        
        
        self.right_layout.addWidget(self.mcs_result_container)
        self.mcs_result_container.show()
        QMessageBox.information(self, "Image Saved", f"The image has been saved to {image_path}")
        QMessageBox.information(self, "Note", "Suppose you try different combinations of scaffolds for MCS. The new image download will replace the old one, so please save it manually if you are satisfied with the results.")


    def close_mcs_result(self):
        if self.mcs_result_container:
            self.right_layout.removeWidget(self.mcs_result_container)
            self.mcs_result_container.deleteLater()
            self.mcs_result_container = None

    
    
    def filter_rule_of_3(self):
        pass
        
    def apply_filter(self, filter_name):
           
        #df = pd.read_csv('fragments.csv')
       
        if filter_name == "Lipinski Rule of 5":
            self.lipin()
        elif filter_name == "Rule of 3 Filter":
            self.ruleof3()
   
    def download_data(self):
        pass
       
    def createRightFrame(self):
        self.right_frame = QWidget()
        self.right_layout = QVBoxLayout(self.right_frame)
        self.table_scroll = QScrollArea()
        self.table_scroll.setWidgetResizable(True)
        self.table_widget = QTableWidget()
        #self.table_widget.setColumnCount(10)
        #self.table_widget.setHorizontalHeaderLabels(["Select","Fragment 1", "Generated Fragments", "Structure"])
        self.table_scroll.setWidget(self.table_widget)
        self.right_layout.addWidget(self.table_scroll)
 


    def ruleof3(self):        
        df = pd.read_csv('fragments.csv')
        results = {"Rule of 3 Filter": []}
        for mol_str in df['Generated Fragments']:
            mol = Chem.MolFromSmiles(mol_str)
            #print(mol)
            if mol:
                if self.check_ruleof3(mol):
                    results["Rule of 3 Filter"].append(mol_str)
        #print("Fragments satisfying Lipinski Rule of 5:")
        filterfragments=[]
       
        for mol_str in results["Rule of 3 Filter"]:
            filterfragments.append(mol_str)
                       
        fl = list(filterfragments)
        unique_fragment_names = self.generate_unique_names(fl)
        self.table_widget.setRowCount(0)    
        for idx, (name, fragment) in enumerate(unique_fragment_names):
            rowPosition = self.table_widget.rowCount()
            self.table_widget.insertRow(rowPosition)
            #-----
            mol_fragment = Chem.MolFromSmiles(fragment)  # Convert fragment SMILES to RDKit molecule
                       
            feature_counts = {
                'Aromatic Rings': Descriptors.NumAromaticRings(mol_fragment),
                'Hydrophobic Atoms': sum(1 for atom in mol_fragment.GetAtoms() if atom.GetAtomicNum() in [6, 1]),
                'Hydrogen Bond Donors': Lipinski.NumHDonors(mol_fragment),
                'Hydrogen Bond Acceptors': Lipinski.NumHAcceptors(mol_fragment),
                'logP_value' : Descriptors.MolLogP(mol_fragment)
                
            }
                   
            checkbox = QCheckBox()
            checkbox.setChecked(False)  # Initially unchecked
            checkbox.fragment = fragment  # Assign fragment to checkbox
            checkbox.stateChanged.connect(self.checkboxStateChanged)  # Connect stateChanged signal to slot function
           
            self.table_widget.setCellWidget(rowPosition, 0, checkbox)
           
            self.table_widget.setItem(rowPosition, 1, QTableWidgetItem(name))
            self.table_widget.setItem(rowPosition, 2, QTableWidgetItem(fragment))
           
            for column, (feature_name, feature_value) in enumerate(feature_counts.items(), start=4):
                self.table_widget.setItem(rowPosition, column, QTableWidgetItem(str(feature_value)))
           
            # Get SMILES structure for the fragment
            smiles_structure = Chem.MolToSmiles(Chem.MolFromSmiles(fragment), isomericSmiles=True)
            self.table_widget.setRowHeight(rowPosition, 60)    
            # Create a QLabel to display the fragment image
            label = QLabel()
            mol = Chem.MolFromSmiles(fragment)
            if mol:
                mol_image = Draw.MolToImage(mol, size=(60, 60))
                qimage = QImage(mol_image.tobytes(), mol_image.width, mol_image.height, QImage.Format_RGB888)
                pixmap = QPixmap.fromImage(qimage)
                label.setPixmap(pixmap)
                self.table_widget.setCellWidget(rowPosition, 3, label)  
                # Set the widget in the cell
           
        self.table_widget.setColumnWidth(0, 40)  
        self.table_widget.setColumnWidth(1, 80)
        self.table_widget.setColumnWidth(2, 200)
        self.table_widget.setColumnWidth(3, 100)
        self.table_widget.setColumnWidth(4, 70)
        self.table_widget.setColumnWidth(5, 70)
        self.table_widget.setColumnWidth(6, 70)
        self.table_widget.setColumnWidth(7, 70)
        self.table_widget.setColumnWidth(8, 70)
        self.table_widget.setColumnWidth(9, 70)
        self.table_widget.setMinimumWidth(900)  
        self.exportTableDataToCSV()    
        self.ruleofthree_fragments = set(unique_fragment_names)
    
    def check_ruleof3(self,mol):
        molecular_weight = Descriptors.ExactMolWt(mol)
        logp = Descriptors.MolLogP(mol)
        h_bond_donor = Descriptors.NumHDonors(mol)
        h_bond_acceptors = Descriptors.NumHAcceptors(mol)
        rotatable_bonds = Descriptors.NumRotatableBonds(mol)
       
        return (molecular_weight <= 300
                and logp <= 3  and
                h_bond_donor <= 3 and
                h_bond_acceptors <= 3 and
                rotatable_bonds <= 3)
                   
    
    
    def exportTableDataToCSV(self):
     headers = []
     for i in range(self.table_widget.columnCount()):
        headers.append(self.table_widget.horizontalHeaderItem(i).text())

     data = []
     for row in range(self.table_widget.rowCount()):
        row_data = []
        for col in range(self.table_widget.columnCount()):
            item = self.table_widget.item(row, col)
            if item is not None:
                row_data.append(item.text())
            else:
                widget = self.table_widget.cellWidget(row, col)
                if isinstance(widget, QCheckBox):
                    row_data.append("Selected" if widget.isChecked() else "Not Selected")
        data.append(row_data)

     with open("fragments.csv", 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(headers)
        csv_writer.writerows(data)
               
   
   
    def updateMolecule(self):
        # Slot to update the molecule visualization when the SMILES string changes
        smiles = self.smiles_textbox.toPlainText().strip()
        if smiles:
            molecule_image = self.loadMolecule(smiles)
            if molecule_image:
                self.image_display.setPixmap(molecule_image.scaled(self.image_display.size(), Qt.KeepAspectRatio))

   
    def loadMolecule(self, smiles):
        # Function to generate a QPixmap from a SMILES string
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            # Generate a 2D image of the molecule
            img = Draw.MolToImage(mol, size=(300, 300))  # Specify the image size here
            buffer = BytesIO()
            img.save(buffer, format="PNG")
            # Convert image to QPixmap
            qt_pixmap = QPixmap()
            qt_pixmap.loadFromData(buffer.getvalue())
            return qt_pixmap
        else:
            QMessageBox.warning(self, "Error", "Invalid SMILES string.")
            return None
       
    def initUI(self):
        # Initialize the UI and create the layout
        self.createLeftFrame()
        self.createRightFrame()
        #self.createImageFrame()  # Create the image frame

        # Create a splitter to manage the layout of frames
        self.splitter = QSplitter(Qt.Horizontal, self)
        self.splitter.addWidget(self.left_frame)
        self.splitter.addWidget(self.right_frame)
        self.splitter.setSizes([340, 400, 650])  # Set initial sizes for the splitter sections

        # Set the central widget and the main window properties
        self.setCentralWidget(self.splitter)
        self.setGeometry(10, 40, 1300, 650)  # Adjust the size and position of the main window
        self.setWindowTitle('Fragment Explorer Version 1.0')
        
        # Show the popup message
        # Show the popup message after the main window is shown
        QTimer.singleShot(500, self.showPopupMessage)


    def showPopupMessage(self):
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Information)
        msg.setText("If you want to try different molecules, please restart the application to avoid incorrect outputs.")
        msg.setWindowTitle("Information")
        msg.setStandardButtons(QMessageBox.Ok)
        msg.exec_()


    def closeFile(self):
        print("Closing the file")  
       
    def showAboutDialog(self):
        print("About the Tool")  # Add your implementation here
   
    def fragbrics(self):
     if hasattr(self, 'smiles_list'):
        smiles_list = self.smiles_list
     else:
        # If not, get SMILES from the text box
        smiles = self.smiles_textbox.toPlainText().strip()
        smiles_list = [smiles]

     for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            fragments_set = BRICS.BRICSDecompose(mol)
            unique_fragment_names = self.generate_unique_names(fragments_set)

     #molecule = self.smiles_textbox.toPlainText().strip()
     #mol = Chem.MolFromSmiles(molecule)
     #fragments_set = BRICS.BRICSDecompose(mol)
     #unique_fragment_names = self.generate_unique_names(fragments_set)

     # Clear existing data in the QTableWidget
     self.table_widget.setRowCount(0)

     # Column titles with full feature names
     column_titles = ["Select", "Fragment 1", "Generated Fragments", "Structure", 'AR', 'HA', 'HBD', 'HBA', 'logP']
     self.table_widget.setColumnCount(len(column_titles))
     self.table_widget.setHorizontalHeaderLabels(column_titles)
 
     # Populate QTableWidget with fragment information
     for idx, (name, fragment) in enumerate(unique_fragment_names):
        rowPosition = self.table_widget.rowCount()
        self.table_widget.insertRow(rowPosition)

        mol_fragment = Chem.MolFromSmiles(fragment)  # Convert fragment SMILES to RDKit molecule

        # Calculate logP and other features
        logP_value = Descriptors.MolLogP(mol_fragment)
        feature_counts = {
            'Aromatic Rings': Descriptors.NumAromaticRings(mol_fragment),
            'Hydrophobic Atoms': sum(1 for atom in mol_fragment.GetAtoms() if atom.GetAtomicNum() in [6, 1]),
            'Hydrogen Bond Donors': Lipinski.NumHDonors(mol_fragment),
            'Hydrogen Bond Acceptors': Lipinski.NumHAcceptors(mol_fragment),
        }

        # Add checkbox for the fragment in the first column
        checkbox = QCheckBox()
        checkbox.setChecked(False)
        checkbox.fragment = fragment
        checkbox.stateChanged.connect(self.checkboxStateChanged)
        self.table_widget.setCellWidget(rowPosition, 0, checkbox)

        # Add fragment name and SMILES
        self.table_widget.setItem(rowPosition, 1, QTableWidgetItem(name))
        self.table_widget.setItem(rowPosition, 2, QTableWidgetItem(fragment))

        # Add feature counts to the table
        for column, (feature_name, feature_value) in enumerate(feature_counts.items(), start=4):
            self.table_widget.setItem(rowPosition, column, QTableWidgetItem(str(feature_value)))

        # Explicitly set the logP value in its column
        self.table_widget.setItem(rowPosition, 8, QTableWidgetItem(str(logP_value)))

        # Get SMILES structure for the fragment
        smiles_structure = Chem.MolToSmiles(Chem.MolFromSmiles(fragment), isomericSmiles=True)
        self.table_widget.setRowHeight(rowPosition, 100)    
        # Create a QLabel to display the fragment image
        label = QLabel()
        mol = Chem.MolFromSmiles(fragment)
        if mol:
                mol_image = Draw.MolToImage(mol, size=(100, 100))
                qimage = QImage(mol_image.tobytes(), mol_image.width, mol_image.height, QImage.Format_RGB888)
                pixmap = QPixmap.fromImage(qimage)
                label.setPixmap(pixmap)
                self.table_widget.setCellWidget(rowPosition, 3, label)  
                # Set the widget in the cell
           
                self.table_widget.setColumnWidth(0, 100)  

        # Set column widths
        self.table_widget.setColumnWidth(0, 100)
        self.table_widget.setColumnWidth(1, 120)
        self.table_widget.setColumnWidth(2, 200)
        self.table_widget.setColumnWidth(3, 100)
        self.table_widget.setColumnWidth(4, 100)
        self.table_widget.setColumnWidth(5, 100)
        self.table_widget.setColumnWidth(6, 100)
        self.table_widget.setColumnWidth(7, 100)
        self.table_widget.setColumnWidth(8, 100)  # logP column width

     self.table_widget.setMinimumWidth(900)
     self.exportTableDataToCSV()
     self.brics_fragments = set(unique_fragment_names)

    def fragrecap(self):
     #molecule = self.smiles_textbox.toPlainText().strip()
     #mol = Chem.MolFromSmiles(molecule)
     #fragment1 = Recap.RecapDecompose(mol)
     #fragments = fragment1.GetAllChildren()
     #fl = list(fragments)
     #unique_fragment_names = self.generate_unique_names(fl)
       
       # Check if SMILES is already loaded from an SDF file
     if hasattr(self, 'smiles_list'):
        smiles_list = self.smiles_list
     else:
        # If not, get SMILES from the text box
        smiles = self.smiles_textbox.toPlainText().strip()
        smiles_list = [smiles]

     for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            fragment1 = Recap.RecapDecompose(mol)
            fragments = fragment1.GetAllChildren()
            fl = list(fragments)
            unique_fragment_names = self.generate_unique_names(fl)
       
     # Clear existing data in the QTableWidget
     self.table_widget.setRowCount(0)
       
           
     # Column titles with full feature names
     column_titles = ["Select", "Fragment 1", "Generated Fragments", "Structure", 'AR', 'HA', 'HBD', 'HBA', 'logP']
     self.table_widget.setColumnCount(len(column_titles))
     self.table_widget.setHorizontalHeaderLabels(column_titles)
 
     #  Populate QTableWidget with fragment information
     for idx, (name, fragment) in enumerate(unique_fragment_names):
        rowPosition = self.table_widget.rowCount()
        self.table_widget.insertRow(rowPosition)

        mol_fragment = Chem.MolFromSmiles(fragment)  # Convert fragment SMILES to RDKit molecule

        # Calculate logP and other features
        logP_value = Descriptors.MolLogP(mol_fragment)
        feature_counts = {
            'Aromatic Rings': Descriptors.NumAromaticRings(mol_fragment),
            'Hydrophobic Atoms': sum(1 for atom in mol_fragment.GetAtoms() if atom.GetAtomicNum() in [6, 1]),
            'Hydrogen Bond Donors': Lipinski.NumHDonors(mol_fragment),
            'Hydrogen Bond Acceptors': Lipinski.NumHAcceptors(mol_fragment),
        }
                   
        # Add checkbox for the fragment in the first column
        checkbox = QCheckBox()
        checkbox.setChecked(False)
        checkbox.fragment = fragment
        checkbox.stateChanged.connect(self.checkboxStateChanged)
        self.table_widget.setCellWidget(rowPosition, 0, checkbox)

        # Add fragment name and SMILES
        self.table_widget.setItem(rowPosition, 1, QTableWidgetItem(name))
        self.table_widget.setItem(rowPosition, 2, QTableWidgetItem(fragment))

        # Add feature counts to the table
        for column, (feature_name, feature_value) in enumerate(feature_counts.items(), start=4):
            self.table_widget.setItem(rowPosition, column, QTableWidgetItem(str(feature_value)))

        # Explicitly set the logP value in its column
        self.table_widget.setItem(rowPosition, 8, QTableWidgetItem(str(logP_value)))
        # Get SMILES structure for the fragment
        smiles_structure = Chem.MolToSmiles(Chem.MolFromSmiles(fragment), isomericSmiles=True)
        self.table_widget.setRowHeight(rowPosition, 100)    
        # Create a QLabel to display the fragment image
        label = QLabel()
        mol = Chem.MolFromSmiles(fragment)
        if mol:
                mol_image = Draw.MolToImage(mol, size=(100, 100))
                qimage = QImage(mol_image.tobytes(), mol_image.width, mol_image.height, QImage.Format_RGB888)
                pixmap = QPixmap.fromImage(qimage)
                label.setPixmap(pixmap)
                self.table_widget.setCellWidget(rowPosition, 3, label)  
                # Set the widget in the cell
           
                self.table_widget.setColumnWidth(0, 100)  

               
        self.table_widget.setColumnWidth(1, 120)
        self.table_widget.setColumnWidth(2, 200)
        self.table_widget.setColumnWidth(3, 100)
        self.table_widget.setColumnWidth(4, 100)
        self.table_widget.setColumnWidth(5, 100)
        self.table_widget.setColumnWidth(6, 100)
        self.table_widget.setColumnWidth(7, 100)
        self.table_widget.setColumnWidth(8, 100)
        self.table_widget.setColumnWidth(9, 100)
       
        self.table_widget.setMinimumWidth(900)  # Adjust width as needed
        self.exportTableDataToCSV()
       
        self.recap_fragments = set(unique_fragment_names)
   

    def smiles_to_sdf(self,smiles):
     # Convert a SMILES string to an SDF format using RDKit
     mol = Chem.MolFromSmiles(smiles)
     if mol is not None:
        AllChem.Compute2DCoords(mol)
        return Chem.MolToMolBlock(mol)
     else:
        return None
   
    def download_sdf(self, smiles, file_name):
     # Remove invalid characters from the file name
     sanitized_file_name = "".join(c for c in file_name if c.isalnum() or c in ['_', '.'])
     # Convert SMILES to SDF and save to a file
     sdf_content = self.smiles_to_sdf(smiles)
     if sdf_content:
        print(f"Saving SDF content to file: {sanitized_file_name}")  # Debug statement
        with open(sanitized_file_name, 'w') as f:
            f.write(sdf_content)
        print("File saved successfully!")  # Debug statement
        QMessageBox.information(self, "File successfully saved", "Verify the working directory you are currently in")

        # Implement the logic to download the file to the user's local system
        # This could involve sending the file to the frontend for download
     else:
        print("Error: Failed to generate SDF content.")  # Debug statement
        QMessageBox.information(self, "File failed to saved", "Verify the working directory you are currently in")
   
   
    def add_sdf_buttons(self, unique_smiles, algorithm, start_row):
     # Add download buttons for SDF files to the table
     for idx, smiles in enumerate(unique_smiles):
        row = start_row + idx
        button = QPushButton("Download")
        button.clicked.connect(lambda _, s=smiles: self.download_sdf(s, f"{s}.sdf"))
        self.unique_table.setCellWidget(row, 3, button)
   
    def compar(self):
     # Extract just the fragment names and SMILES from the stored tuples
     brics_fragments = {frag[1]: frag[0] for frag in self.brics_fragments}  # {SMILES: Name}
     recap_fragments = {frag[1]: frag[0] for frag in self.recap_fragments}  # {SMILES: Name}
    
     unique_brics_smiles = brics_fragments.keys() - recap_fragments.keys()
     unique_recap_smiles = recap_fragments.keys() - brics_fragments.keys()

     # Check if the unique_table already exists in the layout
     if hasattr(self, 'unique_table_container') and self.unique_table_container is not None:
        # Clear the layout that contains the unique_table
        self.right_layout.removeWidget(self.unique_table_container)
        self.unique_table_container.deleteLater()
        self.unique_table_container = None

     # Create a container widget for the unique_table and close button
     self.unique_table_container = QWidget()
     unique_table_layout = QVBoxLayout(self.unique_table_container)

     # Create the close button with an 'X' mark
     close_button = QPushButton('X')
     close_button.setFixedSize(20, 20)  # Set a fixed size for the button (width, height)
     close_button.setStyleSheet("QPushButton { font-size: 10px; }")  # Optional: Adjust font size
     close_button.clicked.connect(self.close_unique_table)

     # Add the close button to the layout
     unique_table_layout.addWidget(close_button)

     # Create the unique_table for displaying unique fragments
     self.unique_table = QTableWidget()
     self.unique_table.setColumnCount(4)
     self.unique_table.setHorizontalHeaderLabels(["Name", "SMILES", "Algorithm", "SDF"])
     self.unique_table.setRowCount(len(unique_brics_smiles) + len(unique_recap_smiles))

     # Populate the unique_table with unique fragments from BRICS
     for idx, smiles in enumerate(unique_brics_smiles):
        name = brics_fragments[smiles]
        self.unique_table.setItem(idx, 0, QTableWidgetItem(name))
        self.unique_table.setItem(idx, 1, QTableWidgetItem(smiles))
        self.unique_table.setItem(idx, 2, QTableWidgetItem("BRICS"))
        # Add download button for each unique fragment
        button = QPushButton("Download")
        button.clicked.connect(lambda _, s=smiles: self.download_sdf(s, f"{name}.sdf"))
        self.unique_table.setCellWidget(idx, 3, button)

     # Populate the unique_table with unique fragments from RECAP
     offset = len(unique_brics_smiles)
     for idx, smiles in enumerate(unique_recap_smiles):
        name = recap_fragments[smiles]
        row = idx + offset
        self.unique_table.setItem(row, 0, QTableWidgetItem(name))
        self.unique_table.setItem(row, 1, QTableWidgetItem(smiles))
        self.unique_table.setItem(row, 2, QTableWidgetItem("RECAP"))
        # Add download button for each unique fragment
        button = QPushButton("Download")
        button.clicked.connect(lambda _, s=smiles: self.download_sdf(s, f"{name}.sdf"))
        self.unique_table.setCellWidget(row, 3, button)

     # Add a label and the unique_table to the unique_table_layout
     unique_table_layout.addWidget(QLabel("Unique Fragments by Algorithm"))
     unique_table_layout.addWidget(self.unique_table)

     # Add the container widget to the right layout
     self.right_layout.addWidget(self.unique_table_container)

    def close_unique_table(self):
     # Slot to handle the close button click
     if self.unique_table_container is not None:
        self.right_layout.removeWidget(self.unique_table_container)
        self.unique_table_container.deleteLater()
        self.unique_table_container = None

    def checkboxStateChanged(self, state):
        checkbox = self.sender()  # Get the checkbox that triggered the signal
        if state == Qt.Checked:
            # Get the associated fragment
            fragment = checkbox.fragment
            # Display the fragment image
            self.displayFragmentImage(fragment)
            self.selected_fragment = checkbox.fragment
            # Uncheck all other checkboxes
            for row in range(self.table_widget.rowCount()):
                if self.table_widget.cellWidget(row, 0) != checkbox:
                    self.table_widget.cellWidget(row, 0).setChecked(False)
                   
                   
     
    def calculate_and_update_features(self, fragment, row):
        # Define your calculation and update logic here
        pass      
   
   
    def visualize_fragments(self, fragments, size=(100, 100), kekulize=True):
        images = []
        for fragment in fragments:
            mol = Chem.MolFromSmiles(fragment)
            images.append(Draw.MolToImage(mol, size=size, kekulize=kekulize))
        return images

    def generate_unique_names(self, fragments):
        unique_names = []
        for idx, fragment in enumerate(fragments):
            unique_name = f"Fragment_{idx}"
            unique_names.append((unique_name, fragment))
        return unique_names

   
    def write_fragments_to_sdf(self, fragments, output_dir='output/'):
        # Create the output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)

        for idx, fragment in enumerate(fragments):
            mol_fragment = Chem.MolFromSmiles(fragment)
            if mol_fragment is not None:
                Chem.MolToMolFile(mol_fragment, os.path.join(output_dir, f"Fragment_{idx}.sdf"))
   
    def displayFragmentImage(self, fragment):
        # Generate a 2D image of the fragment and display it
        mol = Chem.MolFromSmiles(fragment)
        if mol:
            img = Draw.MolToImage(mol, size=(300, 300))  # Specify the image size here
            buffer = BytesIO()
            img.save(buffer, format="PNG")
            # Convert image to QPixmap and display it
            qt_pixmap = QPixmap()
            qt_pixmap.loadFromData(buffer.getvalue())
            self.image_display.setPixmap(qt_pixmap.scaled(self.image_display.size(), Qt.KeepAspectRatio))

    def generateFragments(self):
        print("Generating Fragments")
        # Fetch selected options
        brics_selected = self.left_frame.findChild(QCheckBox, 'BRICS Algorithm').isChecked()
        recap_selected = self.left_frame.findChild(QCheckBox, 'Recap Algorithm').isChecked()

        # Display result in the text area
        if brics_selected and recap_selected:
            self.displayDataInTextArea("Both BRICS Algorithm and Recap Algorithm selected. Hai!")
        elif brics_selected:
            self.displayDataInTextArea("BRICS Algorithm selected. Hai!")
        elif recap_selected:
            self.displayDataInTextArea("Recap Algorithm selected. Hai!")
        else:
            self.displayDataInTextArea("No algorithm selected. Hai!")

    def performStatisticalAnalysis(self):
        print("Performing Statistical Analysis")  # Add your implementation here

    def showGenePlots(self):
        pass  # Add your implementation here
       
    def displayDataInTextArea(self, data):
        self.text_edit.setPlainText(data)  # Set the plain text content of the text area

    def showFileDialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.ReadOnly
        file_dialog = QFileDialog()
        file_dialog.setOptions(options)
        file_dialog.setNameFilter("SDF files  (*.sdf)")

        if file_dialog.exec_() == QFileDialog.Accepted:
            selected_files = file_dialog.selectedFiles()
            #print("Selected file:", selected_files[0])
            try:
                with open(selected_files[0], 'r') as file:
                    sdf_data = file.read()
                # You can now use 'sdf_data' for further processing or display
                print("Ligand Data loaded successfully:\n", sdf_data[:100])  # Displaying only the first 100 characters for illustration
                self.displayDataInTextArea(sdf_data)
            except Exception as e:
                print("Error loading Ligand SDF data:", e)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    tool = Tfrag()
    tool.show()
    sys.exit(app.exec_())
   
    window = MyWindow()
    window.show()
    app.exec_()




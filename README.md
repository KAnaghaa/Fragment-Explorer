**Background**
Drug discovery is a complex and resource-intensive process, often involving the synthesis and screening of large compound libraries.      
**Challenge**
Traditional methods are time-consuming and expensive, necessitating the exploration of more efficient approaches.
**FBDD- Fragment based drug design**
FBDD Offers a promising alternative by breaking down molecules into smaller fragments & systematically building them up. 
**Problem**
Drug discovery is a complex and resource-intensive process, often involving the synthesis and screening of large compound libraries.      

A computational tool called Fragment Explorer was developed to improve and expedite fragment-based drug design procedures. Its primary objective is to identify efficient fragments with the potential to become drug candidates. To achieve this, Fragment Explorer offers a range of powerful features designed to generate, filter, scaffold hopping, maximum common substructure and visualize molecular fragments effectively.

**Principal objective**
Principal Objective: Aims to identify efficient fragments that have the potential to become drug candidates.

**Try it know!**
**Download Fragmeent Explorer 0.1** 
Link: https://github.com/KAnaghaa/Fragment-Explorer/releases/tag/FragmentExplorer
![image](https://github.com/KAnaghaa/Fragment-Explorer/assets/137085789/27e3d09c-9fea-446d-92fb-ecd8362c877f)

**Secondary Objective**
Secondary Objectives:

1.**Fragmentation using BRICS and RECAP Algorithms**:
  a.Uses the BRICS and RECAP algorithms to efficiently fragment molecules. Facilitates the production of various sets of molecular fragments for investigation.
  b.Utilizing advanced algorithms to break down larger molecules into smaller, manageable fragments.
  
2.**Fragment Filtering**:
  a.Creates criteria-based filters, such as the rule of three. Sort fragments according to critical physicochemical and pharmacological characteristics.
  
3.**Scaffold hopping**:
  a.ScaffoldGraph is an open-source cheminformatics library, built using RDKit and NetworkX(study of the structure, dynamics, and functions of complex networks.), for the generation and analysis of scaffold networks and scaffold trees.
  b.Fragments are generated by simplifying a molecule to its core scaffold while removing the side chains and substituents. Its core scaffold while removing side chains and substituents.
  
4.**Maximum common substructure**:
  a.Using the maximum common substructure algorithm from RDKit user can identify common substructures from the scaffold.
  
5.**Visualization Tools**:
  a.Offers 2D representations of fragments to help comprehend their characteristics and relationships. Facilitates better decision-making and data interpretation for researchers.

  ![Screenshot 2024-05-24 120929](https://github.com/KAnaghaa/Fragment-Explorer/assets/137085789/02f4b25d-a247-44c2-872c-e37270fef40f)





**Reference:**
1. Degen, J., Wegscheid-Gerlach, C., Zaliani, A., & Rarey, M. (2008, October 10). On the Art of Compiling and Using “Drug‐Like” Chemical Fragment Spaces. ChemMedChem. https://doi.org/10.1002/cmdc.200800178

2. Lewell, X. Q., Judd, D. B., Watson, S. P., & Hann, M. M. (1998, April 11). RECAPRetrosynthetic Combinatorial Analysis Procedure: A Powerful New Technique for Identifying Privileged Molecular Fragments with Useful Applications in Combinatorial Chemistry. Journal of Chemical Information and Computer Sciences. https://doi.org/10.1021/ci970429i

3.Cao, Y., Jiang, T., & Girke, T. (2008, July 1). A maximum common substructure-based algorithm for searching and predicting drug-like compounds. Bioinformatics. https://doi.org/10.1093/bioinformatics/btn186

4.Kralj, S., Jukič, M., & Bren, U. (2023, April 18). Molecular Filters in Medicinal Chemistry. Encyclopedia. https://doi.org/10.3390/encyclopedia3020035

5. Li, Q. (2020, August 5). Application of Fragment-Based Drug Discovery to Versatile Targets. Frontiers in Molecular Biosciences. https://doi.org/10.3389/fmolb.2020.00180

6.Ivanov, N. N., Shulga, D. A., & Palyulin, V. A. (2023, May 24). Decomposition of Small Molecules for Fragment-Based Drug Design. Biophysica. https://doi.org/10.3390/biophysica3020024

7.Konteatis, Z. (2021, March 26). What makes a good fragment in fragment-based drug discovery? Expert Opinion on Drug Discovery. https://doi.org/10.1080/17460441.2021.1905629

8.Scott, O. B., & Chan, A. W. E. (2020, March 31). ScaffoldGraph: an open-source library for the generation and analysis of molecular scaffold networks and scaffold trees. Bioinformatics. https://doi.org/10.1093/bioinformatics/btaa219

**Acknowledgement**

1.**SASTRA DEEMED UNIVERSITY**: I would like to express my sincere gratitude to **Mr. Udayakumar.M (Asst. Professor - III, School of Chemical and Biotechnology)** for his invaluable guidance throughout the development of Fragment Explorer: A Comprehensive Tool for Fragment-Based Drug Design - PyQt Application. His insights and support were instrumental in shaping this project into a successful endeavor.

I would like to thank **SASTRA DEEMED UNIVERSITY** for its infrastructure and facilities.

2.**ScaffoldGraph from the paper**:  Scott, O. B., & Chan, A. W. E. (2020, March 31). ScaffoldGraph: an open-source library for the generation and analysis of molecular scaffold networks and scaffold trees. Bioinformatics. https://doi.org/10.1093/bioinformatics/btaa219
   link: https://github.com/UCLCheminformatics/ScaffoldGraph

3.**RDKit**: Open-source cheminformatics. https://www.rdkit.org/






# 🦟 Estimating Recent Effective Population size from Resistance Loci in _Anopheles gambiae_

We present MEOWSS — _Multi-origin Estimation Of Widespread Soft Sweeps_ <br/>
A simple heuristic method that applies time to the most recent common ancestor (TMRCA)-based hierarchical clustering on haplotypes sampled from populations undergoing soft selective sweeps.
It is designed to estimate the number of independent origins (𝜂) and infer very recent effective population sizes (N<sub>e</sub>) using the semi-deterministic theory [(Khatri & Burt, 2019)](https://academic.oup.com/mbe/article/36/9/2040/5435963) <br/>

🐈‍⬛ Explanation for why MEOWSS is called MEOWSS... 🐈‍⬛<br/>
- _Multi-origin_ ➡ refers to the fact that soft sweeps arise when a beneficial mutation appears independently on multiple haplotypes
- _Estimation_ ➡  indicates that the method provides a quantitative heuristic to infer the number of origins
- _of Widespread Soft Sweeps_ ➡ targets soft selective sweeps!
- its a super cute name _duuh_...

## 📖 About

- **University:** Imperial College London
- **Course:** BSc Biochemistry 
- **Supervisor:** Dr Bhavin Khatri

## 🪜 MEOWSS Workflow
a brief outline of the experimental workflow, more details on [my report](https://docs.google.com/document/d/151sB2UNQ25oAq-I8lZecM9KjynHKLwHxQE8o99DtFuk/edit?usp=sharing). Most scripts were run on the Imperial High Performance Cluster (HPC) using a `.sh` file.
### 1️⃣ Soft Sweep Simulations in [SLiM 4](https://messerlab.org/slim/)
1. Create burn-in populations (store in `.trees` files)
2. Start Soft Sweep simulation with the burn-ins<br/>
    📂 [SS_all.slim](https://github.com/kitku15/MEOWSS/blob/master/SLIM_SIMULATION/SS_all.slim)<br/>
    📂 [run_simulations.sh](https://github.com/kitku15/MEOWSS/blob/master/HPC_JOB_FILES/run_simulations.sh)
3. outputs: `.trees`, `.vcf`, `.sample_ID.txt`, `log.csv`<br/>
    📂 [TREES](https://github.com/kitku15/MEOWSS/tree/master/EXAMPLE_FILES/TREES)<br/>
    📂 [VCF_NEW](https://github.com/kitku15/MEOWSS/tree/master/EXAMPLE_FILES/VCF_NEW) _this is post preprocessing btw!_<br/>
    📂 [sample_ID](https://github.com/kitku15/MEOWSS/tree/master/EXAMPLE_FILES/sample_ID)

💡*reading [SLiM](https://github.com/MesserLab/SLiM/releases/download/v5.0/SLiM_Manual.pdf) and [Eidos](https://github.com/MesserLab/SLiM/releases/download/v5.0/Eidos_Manual.pdf) manuals helped! They also have easy to follow along video tutorials on their [website](https://messerlab.org/slim/#Workshops)!*

### 2️⃣ VCF Preprocessing
1. Remove all m0 mutations at the sweep site <br/>
    📂 [fix_vcf.py](https://github.com/kitku15/MEOWSS/blob/master/SLIM_SIMULATION/fix_vcf.py)
2. convert .vcf into .vcz <br/>
    📂 [run_vcf_to_zarr.sh](https://github.com/kitku15/MEOWSS/blob/master/HPC_JOB_FILES/run_vcf_to_zarr.sh)
### 3️⃣ MEOWSS TMRCA-based Haplotype Clustering
1. input: `.vcz`, threshold for t<sub>hom</sub> & smoothing points<br/>
    📂 [get_active_zars.py](https://github.com/kitku15/MEOWSS/blob/master/tmrca/get_active_zars.py)<br/>
    📂 [active_zars.txt](https://github.com/kitku15/MEOWSS/blob/master/tmrca/active_zars.txt)
2. TMRCA Clustering <br/>
    📂 [tmrca_MEOWSS.py](https://github.com/kitku15/MEOWSS/blob/master/tmrca/tmrca_MEOWSS.py)<br/>
    📂 [run_tmrca_MEOWSS.sh](https://github.com/kitku15/MEOWSS/blob/master/HPC_JOB_FILES/run_tmrca_MEOWSS.sh)
3. output: `Z.npy`, `cols.npy`, unlabelled dendrogram<br/>
    📂 [Z.npy](https://github.com/kitku15/MEOWSS/blob/master/EXAMPLE_FILES/TMRCA/0.5/100/1_13_1000_0.00025_2.5e-05_0.5_100_0.97_10_100_Z.npy)<br/>
    📂 [cols.npy](https://github.com/kitku15/MEOWSS/blob/master/EXAMPLE_FILES/TMRCA/0.5/100/1_13_1000_0.00025_2.5e-05_0.5_100_0.97_10_100_cols.npy)<br/>
    📂 [unlabelled dendrogram](https://github.com/kitku15/MEOWSS/blob/master/EXAMPLE_FILES/TMRCA/0.5/100/1_13_1000_0.00025_2.5e-05_0.5_100_0.97_10_100.pdf)<br/>
    
💡*reading methods section of my report will help you understand how it works*
### 4️⃣ Analysis
1. input: `Z.npy`, `cols.npy`, `.trees` from above<br/>
2. scores dendrograms, labels dendrograms by true haplotype origin, Determining 𝜂 for t<sub>clus</sub> = 1-4<br/>
    📂 [analysis.py](https://github.com/kitku15/MEOWSS/blob/master/tmrca/analysis.py)<br/>
    📂 [analysis_trees.py](https://github.com/kitku15/MEOWSS/blob/master/tmrca/analysis_trees.py)<br/>
    output to meta files<br/>
    📂 [METAS](https://github.com/kitku15/MEOWSS/tree/master/METAS)
5. Make perfect tree linkages<br/>
    📂 [plot_my_trees](https://github.com/kitku15/MEOWSS/tree/master/plot_my_trees)<br/>
    📂 [run_perfectTREES.sh](https://github.com/kitku15/MEOWSS/blob/master/HPC_JOB_FILES/run_perfectTREES.sh)
5. get cophenetic correlation coefficient between Z linkage from true genealogy and MEOWSS-inferred genealogy<br/>
    📂 [Z_analysis.py](https://github.com/kitku15/MEOWSS/blob/master/tmrca/Z_analysis.py)<br/>
    📂 [run_Z_analysis.sh](https://github.com/kitku15/MEOWSS/blob/master/HPC_JOB_FILES/run_Z_analysis.sh)<br/>
6. Calculate N<sub>e</sub> by maximum likelihood estimation <br/>
    📂 [MLE_plot.py](https://github.com/kitku15/MEOWSS/blob/master/tmrca/MLE_plot.py)<br/>
    📂 [MLE_plots](https://github.com/kitku15/MEOWSS/tree/master/tmrca/MLE_plots)<br/>

##  🎶 SINGER
Simulated VCF files were processed with SINGER (Deng et al., 2024), producing an ARG in .trees format (tskit). For each simulation, the tree at the sweep site was used to construct a linkage matrix, which was then used to generate a dendrogram as described above.<br/>
📂[get_vcfs_for_SINGER.txt](https://github.com/kitku15/MEOWSS/blob/master/tmrca/get_vcfs_for_SINGER.txt)<br/>
📂[run_singer.sh](https://github.com/kitku15/MEOWSS/blob/master/HPC_JOB_FILES/run_singer.sh)<br/>
📂[SINGER](https://github.com/kitku15/MEOWSS/tree/master/EXAMPLE_FILES/SINGER)<br/>
📂[singer_meta.csv](https://github.com/kitku15/MEOWSS/blob/master/METAS/singer_meta.csv)<br/>

##  🧬 Phasedibd
Simulated VCF files were processed with PhasedIBD (Freyman et al., 2020) to get IBD segments around the sweep site which were used as the SHLs in a similar tmrca clustering method as MEOWSS.<br/>
📂[tmrca_phaseibd.py](https://github.com/kitku15/MEOWSS/blob/master/tmrca/tmrca_phaseibd.py)<br/>
📂[PHASEIBD](https://github.com/kitku15/MEOWSS/tree/master/EXAMPLE_FILES/PHASEIBD)<br/>
📂[phasedibd_meta.csv](https://github.com/kitku15/MEOWSS/blob/master/METAS/phasedibd_meta.csv)<br/>

## 🏅 Acknowledgements
To Dr Bhavin Khatri and Josh Reynolds for their supervision. I acknowledge Cheyanne Seah, Anushka Thawani and Theo Hemmant for earlier contributions to the Python code implemented in MEOWSS and experimental protocol, and the Imperial College Research Computing Service for computational support.

## 🫂 To the next Student!

Inheriting this project was super tricky but it gets better (and more fun!) as you progress! just keep going and dont give up! :3 Spend the first month understanding the research area itself and the scripts. I've put down future directions in the discussions section of the report as potential ideas for you to continue!

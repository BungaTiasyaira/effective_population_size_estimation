# 🦟 Estimating Recent Effective Population size from Resistance Loci in _Anopheles gambiae_

We present MEOWSS — _Multi-origin Estimation Of Widespread Soft Sweeps_ <br/>
A simple heuristic method that applies time to the most recent common ancestor (TMRCA)-based hierarchical clustering on haplotypes sampled from populations undergoing soft selective sweeps.
It is designed to estimate the number of independent origins (𝜂) and infer very recent effective population sizes (N<sub>e</sub>) using the semi-deterministic theory [(Khatri & Burt, 2019)](https://academic.oup.com/mbe/article/36/9/2040/5435963) <br/>

## 📖 About

- **University:** Imperial College London
- **Course:** BSc Biochemistry 
- **Supervisor:** Dr Bhavin Khatri

## 🪜 Workflow
a brief outline of the experimental workflow, more details on [my report](https://docs.google.com/document/d/151sB2UNQ25oAq-I8lZecM9KjynHKLwHxQE8o99DtFuk/edit?usp=sharing)
### 1️⃣ Soft Sweep Simulations in [SLiM 4](https://messerlab.org/slim/)
1. Create burn-in populations (store in `.trees` files)
2. Start Soft Sweep simulation with the burn-ins
3. outputs: `.trees`, `.vcf`, `.sample_ID.txt`
### 2️⃣ VCF Preprocessing
1. Remove all m0 mutations at the sweep site
2. convert .vcf into .vcz
### 3️⃣ MEOWSS TMRCA-based Haplotype Clustering
1. input: `.vcz`, threshold for t<sub>hom</sub> & smoothing points
2. Estimating Shared Haplotype Length 
3. Calculating TMRCA of Haplotype Pairs
4. Haplotype Clustering 
5. output: `Z.npy`, `cols.npy`, dendrogram
### 4️⃣ Analysis
1. intput: `Z.npy`, `cols.npy`, `.trees`
2. scores dendrograms
3. labels dendrograms by true haplotype origin
4. Determining 𝜂 for t<sub>clus</sub> = 1-4
5. get cophenetic correlation coefficient between Z linkage from true genealogy and MEOWSS-inferred genealogy
6. Calculate N<sub>e</sub> by maximum likelihood estimation 


## 🫂 To the next Student!

Inheriting this project was super tricky but it gets better (and more fun!) as you progress! just keep going and dont give up! :3

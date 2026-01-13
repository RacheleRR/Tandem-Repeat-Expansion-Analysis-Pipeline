

## 🚀 Installation
### 1.Setup
```bash
# 1. Clone repository
git clone https://github.com/yourusername/ehdn-pipeline.git
cd ehdn-pipeline

# 2. Run automated setup (downloads tools, resources, helper scripts)
bash setup/setup.sh

```
### 2. Manuali download annovar files 
register and download annovar from this website 
https://www.openbioinformatics.org/annovar/annovar_download_form.php
or this https://annovar.openbioinformatics.org/en/latest/user-guide/download/
```bash
tar xvfz annovar.latest.tar.gz 
mv annovar helper/annovar 
```

### 3. Download gmt from gprofiler 
If you desire to do the network analysis you need to go onto https://biit.cs.ut.ee/gprofiler/gost and to Data source select the data sources you would like to use for your network analysis and download the combined names gmt 

```bash
mv combined_names.gmt  resources/ 
```


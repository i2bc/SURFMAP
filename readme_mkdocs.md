# MkDocs

### Official documentation link

https://www.mkdocs.org/getting-started/


## How it works

[Writing your docs](https://www.mkdocs.org/user-guide/writing-your-docs/)

Based on:
 -`mkdocs.yml`: high-level configuration file structuring all your documentation pages.
 - `docs/`: main directory the mkdocs.yml file automatically refers to. It contains all the documentation pages and other accessory files/directory to customize the pages (css, javascript) or add images.
  
 
###  Example of docs/
<pre><font color="#3465A4"><b>docs/</b></font>
├── <font color="#3465A4"><b>css/</b></font>
├── <font color="#4E9A06"><b>download.md</b></font>
├── <font color="#4E9A06"><b>How_it_works_orfold.md</b></font>
├── <font color="#3465A4"><b>img/</b></font>
├── <font color="#4E9A06"><b>index.md</b></font>
├── <font color="#4E9A06"><b>index_orfmap.md</b></font>
├── <font color="#4E9A06"><b>index_orfmine.md</b></font>
├── <font color="#4E9A06"><b>index_ori.md</b></font>
├── <font color="#4E9A06"><b>installation.md</b></font>
├── <font color="#3465A4"><b>js/</b></font>
├── <font color="#4E9A06"><b>Objective_orfold.md</b></font>
├── <font color="#4E9A06"><b>orfget_param.md</b></font>
├── <font color="#4E9A06"><b>orfget_run.md</b></font>
├── <font color="#4E9A06"><b>orfmine_citations.md</b></font>
├── <font color="#4E9A06"><b>orfmine_installation.md</b></font>
├── <font color="#4E9A06"><b>orfmine_overview.md</b></font>
├── <font color="#4E9A06"><b>orfmine_quickstart.md</b></font>
├── <font color="#4E9A06"><b>orftrack_annotation.md</b></font>
├── <font color="#4E9A06"><b>orftrack_description.md</b></font>
├── <font color="#4E9A06"><b>orftrack_input.md</b></font>
├── <font color="#4E9A06"><b>orftrack_orfdef.md</b></font>
├── <font color="#4E9A06"><b>orftrack_orf_extraction.md</b></font>
├── <font color="#4E9A06"><b>orftrack_overlap.md</b></font>
├── <font color="#4E9A06"><b>orftrack_param.md</b></font>
├── <font color="#4E9A06"><b>orftrack_run.md</b></font>
├── <font color="#4E9A06"><b>parameters_orfold.md</b></font>
├── <font color="#4E9A06"><b>parameters_orfplot.md</b></font>
├── <font color="#4E9A06"><b>Plot_orfold.md</b></font>
├── <font color="#4E9A06"><b>Run_orfold_advanced.md</b></font>
├── <font color="#4E9A06"><b>Run_orfold.md</b></font>
└── <font color="#4E9A06"><b>virtualenv.md</b></font>
</pre>

 
###  Example of mkdocs.yml

Please note that it's a YAML file, so indentation is important.


```
site_name: ORFmine Documentation
nav:
    - Home: index.md
    - ORFmine:
      - Installation: orfmine_installation.md
      - Quick start: orfmine_quickstart.md
      - Citations: orfmine_citations.md
    - ORFtrack:
      - ORFtrack presentation:
        - Aims and general description: orftrack_description.md
        - ORF definition: orftrack_orfdef.md
        - ORF annotation: orftrack_annotation.md
        - Overlap definition: orftrack_overlap.md
        - ORF extraction: orftrack_orf_extraction.md
      - Usage:
        - Inputs: orftrack_input.md
        - ORF annotation with ORFtrack: orftrack_run.md
        - ORF extraction with ORFget: orfget_run.md
      - Parameters:
        - ORFtrack parameters: orftrack_param.md
        - ORFget parameters: orfget_param.md
    - ORFold:
      - ORFold presentation:
        - Objective: Objective_orfold.md
        - How it works?: How_it_works_orfold.md
      - Usage:
        - Basic run ORFold: Run_orfold.md
        - Plots with ORFplot: Plot_orfold.md
        - Advanced run ORFold: Run_orfold_advanced.md
      - Parameters:
        - ORFold parameters: parameters_orfold.md
        - ORFplot parameters: parameters_orfplot.md
theme: readthedocs
markdown_extensions:
    - footnotes
    - attr_list
    - toc:
        permalink: "#"
use_directory_urls: false
extra_css:
    - css/extra.css
extra_javascript:
    - js/extra.js
docs_dir: docs
```

Above is shown a file with 7 mkdocs setting keys: site_name, nav, theme, markdown_extensions, use_directory_urls, extra_css, extra_javascript.

- `site name:` is followed by the name/title you want to give to your documentation site
- `nav:` lists different navigation link each referring to a markdown page. Those pages are your documentation page (they can be converted in html pages later by mkdocs). For the example above,
if a user click on the `Quick start` link, he will be redirected to the page `orfmine_quickstart.md`. 
- `docs_dir:` indicates the directory containing the documentation source markdown files (and others such as `css/extra.css`). This can either be a relative directory, in which case it is resolved relative to `mkdocs.yml`, or it can be an absolute directory path from the root of your local file system.


Please check the mkdocs official documentation for more information on other [mkdocs setting](https://www.mkdocs.org/user-guide/configuration/) keys that can be used.  


## Preview your documentation in live

MkDocs comes with a built-in dev-server that lets you preview your documentation as you work on it. Make sure you're in the same directory as the mkdocs.yml configuration file, and then start the server by running the command:
```
mkdocs serve
```

<pre>$ mkdocs serve
INFO     -  Building documentation...
INFO     -  Cleaning site directory
INFO     -  The following pages exist in the docs directory, but are not included in the &quot;nav&quot;
            configuration:
              - download.md
              - index_orfmap.md
              - index_orfmine.md
              - index_ori.md
              - installation.md
              - orfmine_overview.md
              - virtualenv.md
INFO     -  Documentation built in 0.17 seconds
INFO     -  [10:49:54] Watching paths for changes: &apos;docs&apos;, &apos;mkdocs.yml&apos;
INFO     -  [10:49:54] Serving on http://127.0.0.1:8000/
</pre>

Open up http://127.0.0.1:8000/ in your browser, and you'll see the default home page being displayed.

The dev-server also supports auto-reloading, and will rebuild your documentation whenever anything in the configuration file, documentation directory, or theme directory changes.





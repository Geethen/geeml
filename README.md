<div id="top"></div>

<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->
[![image](https://img.shields.io/pypi/v/geeml.svg)](https://pypi.python.org/pypi/geeml)
[![image](https://img.shields.io/conda/vn/conda-forge/geeml.svg)](https://anaconda.org/conda-forge/geeml)
[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]


<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/Geethen/geeml">
    <img src="./images/logo2_GEEML.png" alt="Logo" width="400" height="400">
  </a>

<h3 align="center">GEEML: Google Earth Engine Machine learning</h3>

  <p align="center">
    A python package to extract gee data for machine learning.
    <br />
    <a href="https://Geethen.github.io/geeml"><strong>Explore the documentation »</strong></a>
    <br />
    <br />
    <a href="https://github.com/Geethen/geeml">View Demo</a>
    ·
    <a href="https://github.com/Geethen/geeml/issues">Report Bug</a>
    ·
    <a href="https://github.com/Geethen/geeml/issues">Request Feature</a>
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Basic Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project
This python package makes it easier to extract satellite data from Google Earth Engine using parallel processing and the Google Earth Engine high volume end point.

In its current state it supports the extraction of data for traditional machine learning (tabular data) in the form of csv's and the extraction of GeoTiff image patches for Deep Neural Networks.

#### Motivation
The Machine learning capabilities in the GEE JS code editor remain limited. For example, there is no support for XGBoost, LightGBM, NGBoost, etc. Moreover, the python ecosystem has much more support for training, valdation and hyperparameter tuning. However, for this functionality to be leveraged, data needs to be downloaded locally or stored in Google Drive or Google Cloud Storage to benefit from the Machine learning python ecosystem. Therfore, this package aims to make it easier and faster to download GEE-processed data in a machine learning-ready format. 

**Features**
* Parallel export of images or sparse images (for example, GEDI).
* Export raster values at points or polygons (ee.FeatureCollection).
* Summarise raster data within polygons (ee.FeatureCollections).
* Extract both tabular and Deep Neural Network (DNN) type datasets.

<!-- GETTING STARTED -->
## Getting Started

### Installation
To install this package:

1. pip 
   ```sh
   pip install geeml
   ```
2. Build from source (latest version)
   ```sh
   pip install git+https://github.com/Geethen/geeml.git
   ```

<p align="right">(<a href="#top">back to top</a>)</p>


<!-- USAGE EXAMPLES -->
## Basic usage
Download the NASADEM elevation data for Kenya.

   ```python
  #import packages
  import ee
  from geeml.utils import getCountry
  from geeml.extract import extractor

  # Authenticate GEE
  ee.Authenticate()
  # Initialize GEE with high-volume end-point
  ee.Initialize(opt_url='https://earthengine-highvolume.googleapis.com')
   
  # Import datasets from GEE
  nasadem = ee.Image("NASA/NASADEM_HGT/001")
  # A point in Kenya
  poi = ee.Geometry.Point([37.857884,-0.002197])
  kenya = getCountry(poi)

  # Download directory
  dd = '/content/drive/MyDrive/geeml_example'

  # Prepare for data extraction
  trialExtractor = extractor(covariates=nasadem, aoi = kenya, scale= 5000, dd= dd)

  # Extract data
  trialExtractor.extractAoi()
   ```

_For more examples, please refer to the [Documentation](https://geethen.github.io/geeml/notebooks/example/)_

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- ROADMAP -->
## Roadmap

- [ ] Support the export of additional formats (TFrecords)
- [ ] Download data from GEE based on local shapefiles
- [ ] Add more examples for using the package

See the [open issues](https://github.com/Geethen/geeml/issues) for a full list of proposed features (and known issues).

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to learn, inspire, and create. Any contributions you make are **greatly appreciated**.

If you have a suggestion that would make this better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement".
Don't forget to give the project a star! Thanks again!

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- CONTACT -->
## Contact

Geethen Singh - [@Geethen](https://twitter.com/Geethen) - geethen.singh@gmail.com

Project Link: [https://github.com/Geethen/geeml](https://github.com/Geethen/geeml)

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- ACKNOWLEDGMENTS -->
## Acknowledgments

* [Natural State](https://www.naturalstate.org/)
* [University of the Witwatersrand](https://www.wits.ac.za/)
* [Fitz patrick centre for african ornithology](http://www.fitzpatrick.uct.ac.za/)

This package was created with [Cookiecutter](https://github.com/cookiecutter/cookiecutter) and the [giswqs/pypackage](https://github.com/giswqs/pypackage) project template.

This package uses the [geedim](https://pypi.org/project/geedim/) package for extracting image data at an AOI

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/Geethen/geeml.svg?style=for-the-badge
[contributors-url]: https://github.com/Geethen/geeml/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/Geethen/geeml.svg?style=for-the-badge
[forks-url]: https://github.com/Geethen/geeml/network/members
[stars-shield]: https://img.shields.io/github/stars/Geethen/geeml.svg?style=for-the-badge
[stars-url]: https://github.com/Geethen/geeml/stargazers
[issues-shield]: https://img.shields.io/github/issues/Geethen/geeml.svg?style=for-the-badge
[issues-url]: https://github.com/Geethen/geeml/issues
[license-shield]: https://img.shields.io/github/license/Geethen/geeml.svg?style=for-the-badge
[license-url]: https://github.com/Geethen/geeml/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/linkedin_username
[product-screenshot]: images/screenshot.png
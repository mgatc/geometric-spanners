<div id="top"></div>




<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->
[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]



<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/mgatc/geometric-spanners">
<h3 align="center">Bounded-degree Plane Geometric Spanner Testbed</h3>
  </a>

  <p>
    The construction of bounded-degree plane geometric spanners has been a focus of interest  since 2002 when Bose, Gudmundsson, and Smid proposed the first algorithm to construct such spanners. To date, eleven algorithms have been designed with various trade-offs in degree and stretch factor. We have carefully implemented these sophisticated algorithms in \textsf{C}\texttt{++} using the  \textsf{CGAL} library and experimented with them using large synthetic and real-world pointsets. Our  experiments  have revealed their practical behavior on large pointsets and their real-world efficacy. We share the implementations here for broader uses and future research.  
  </p>
  <p>
    We present a simple practical algorithm, named \textsc{AppxStretchFactor}, that can estimate  stretch factors (obtains  lower bounds on the exact stretch factors) of geometric spanners fast -- a challenging problem for which no practical algorithm is known yet. In our experiments with bounded-degree plane geometric spanners, we find that \textsc{AppxStretchFactor} estimates stretch factors almost precisely and gives linear runtime performance for the pointset distributions considered. We also found it to be much faster than the naive Dijkstra-based algorithm for calculating stretch factors. Further, it can be easily parallelized, making it very useful for estimating stretch factors of large spanners.
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->



## Built With

* [GCC](https://gcc.gnu.org/)
* [C++17](https://en.cppreference.com/w/cpp/17)
* [CGAL](https://www.cgal.org/)
* [Boost](https://www.boost.org/)
* [GMP](https://gmplib.org/)
* [MPFR](https://www.mpfr.org/)

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- GETTING STARTED -->
## Getting Started


### Prerequisites

Ensure the dependencies listed above are available. 

Additionally, `texlive-full` is required for visualizations

### Installation

1. Clone the repo
   ```sh
   git clone https://github.com/mgatc/geometric-spanners.git
   ```
2. Attempt to build
   ```sh
   cmake --build ./geometric-spanners/cmake-build-debug --target spanners
   ```
3. If unsuccessful, try configuring ```CMakeLists.txt``` to find the dependencies.

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- USAGE EXAMPLES -->
## Usage

* `./spanners` to view a brief usage summary

### Viewing visualizations

* `./spanners [n]` to produce a visualization of the algorithm configured in `src/Scratch.h` with `n` points
* `./spanners [filename.xy]` to produce a visualization of the algorithm configured in `src/Scratch.h` with the points contained in `filename.xy` (the file extension **must** be `.xy`)


### Running experiments

* `./spanners [filename.xml] [r]` to run an experiment on real-world point sets defined in `filename.xml` with `r` iterations
*  `./spanners [r] [start n] [end n] [increment]` to run an experiment on eight synthetic distributions with `r` iterations

<p align="right">(<a href="#top">back to top</a>)</p>




<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- ACKNOWLEDGMENTS -->
## Acknowledgments

* [David Wisnosky](https://github.com/Wisno33)
* [Lucas Mougeot](https://github.com/lucasfuturist)
* [Anirban Ghosh](https://github.com/ghoshanirban)
* [Fred Anderson](https://github.com/TheDKG)

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/mgatc/geometric-spanners.svg?style=for-the-badge
[contributors-url]: https://github.com/mgatc/geometric-spanners/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/mgatc/geometric-spanners.svg?style=for-the-badge
[forks-url]: https://github.com/mgatc/geometric-spanners/network/members
[stars-shield]: https://img.shields.io/github/stars/mgatc/geometric-spanners.svg?style=for-the-badge
[stars-url]: https://github.com/mgatc/geometric-spanners/stargazers
[issues-shield]: https://img.shields.io/github/issues/mgatc/geometric-spanners.svg?style=for-the-badge
[issues-url]: https://github.com/mgatc/geometric-spanners/issues
[license-shield]: https://img.shields.io/github/license/mgatc/geometric-spanners.svg?style=for-the-badge
[license-url]: https://github.com/mgatc/geometric-spanners/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/matt-graham-368509106
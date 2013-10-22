Diffusion Tensor Imaging (DTI) operator using a "spectral quaternion" decomposition.
====================================================================================

The repository contains Anne Collard's work on "diffusion tensor imaging". The central part is the 'difftensor' matlab object defined using a spectral-quaternion decomposition. This matlab object lets you manipulate diffusion tensors by overloading standard operations such as the distance between tensors, interpolation and averaging. See the object's help for more details. Secondly the repository contains a set of functions illustrating how the difftensor object can be employed to implement practical operations on DTI's.


References
----------
Anne Collard, PhD. thesis,"Geometric statistical processing of brain diffusion tensor images", 2013, http://hdl.handle.net/2268/156140.

Contacts:
--------

- Dr. Anne Collard  
  Department of Electrical Engineering and Computer Science, University of Liège, Belgium  
  anne.collard_at_ulg.ac.be  
- Prof. Christophe Phillips  
  Cyclotron Research Centre, University of Liège, Belgium  
  c.phillips_at_ulg.ac.be
- Prof. Rodolphe Sepulchre  
  Department of Electrical Engineering and Computer Science, University of Liège, Belgium  
  r.sepulchre_at_ulg.ac.be  
  and  
  Department of Engineering, University of Cambridge, UK  
  r.sepulchre_at_eng.cam.ac.uk

Thesis abstract:
---------------
Nowadays, the functioning of the human brain is still one of the biggest mysteries of humanity. The multiple holes in the understanding of the human brain explain why an intensification of brain-oriented research can be observed since a few years.

One of the most recent techniques to better understand the brain is Diffusion Tensor Imaging (DTI), a noninvasive imaging modality that provides information about orientation of nervous fibers, and their spatial density, with a high resolution. The particular nature of DTI images makes them multi-valued. Their processing therefore requires to adapt state-of-the-art techniques, which are fundamentally tailored to scalar-valued images. 

The objective of this PhD thesis is to develop a novel framework for the processing of tensor diffusion images. The focus is threefold: first, we adopt a Riemannian geometric framework to generalize image processing from linear to nonlinear spaces. Second, we aim at developing a processing framework that retains the physical information of measurement data. Thirdly, the proposed algorithm must be computationally efficient in order to scale with the data size of clinical applications. 

The main contribution of this thesis is the development of a novel processing method, which has the particularity to preserve the important features of diffusion tensors, while being computationally affordable. This technique is based on the decoupling between the two types of information conveyed by tensors: the diffusion intensity on one hand, and the orientation of diffusion on the other hand. Moreover, the computational cost is limited thanks to the use of unit quaternions to represent tensors orientation. 

Another contribution of the thesis lies in the development of a statistical method for group comparison. This method uses the notion of similarity measure between the values, a notion that can be defined for multi-valued images, and which enables to reduce the computational cost. The use, for the statistical tests, of the similarity measure associated to our framework turns out to be efficient and informative. 

The study of geometric methods for multi-valued images together with the study of potential applications of diffusion tensor images have enabled the introduction of a novel framework, which is particularly appropriate for those images. The basic operations developed in the thesis open the way to more sophisticated processing algorithms, while ensuring the preservation of the main information associated to the tensors.

DISCLAIMER:
----------
This repository content is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free oftware Foundation,  either version 2 of  the License,  or (at  your option) any later version. 
This code is distributed in the hope  that it will be  useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A  PARTICULAR PURPOSE.  See the  GNU General Public License for more  details.
A copy of the  GNU General Public License is available in this repository. If not, see <http://www.gnu.org/licenses/>.


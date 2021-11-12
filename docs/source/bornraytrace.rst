<!DOCTYPE html>
<html class="writer-html5" lang="en" >
<head>
  <meta charset="utf-8" /><meta name="generator" content="Docutils 0.17.1: http://docutils.sourceforge.net/" />

  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>bornraytrace package &mdash; BornRaytrace  documentation</title>
      <link rel="stylesheet" href="_static/pygments.css" type="text/css" />
      <link rel="stylesheet" href="_static/css/theme.css" type="text/css" />
  <!--[if lt IE 9]>
    <script src="_static/js/html5shiv.min.js"></script>
  <![endif]-->
  
        <script data-url_root="./" id="documentation_options" src="_static/documentation_options.js"></script>
        <script src="_static/jquery.js"></script>
        <script src="_static/underscore.js"></script>
        <script src="_static/doctools.js"></script>
    <script src="_static/js/theme.js"></script>
    <link rel="index" title="Index" href="genindex.html" />
    <link rel="search" title="Search" href="search.html" />
    <link rel="prev" title="bornraytrace" href="modules.html" /> 
</head>

<body class="wy-body-for-nav"> 
  <div class="wy-grid-for-nav">
    <nav data-toggle="wy-nav-shift" class="wy-nav-side">
      <div class="wy-side-scroll">
        <div class="wy-side-nav-search" >
            <a href="index.html" class="icon icon-home"> BornRaytrace
          </a>
<div role="search">
  <form id="rtd-search-form" class="wy-form" action="search.html" method="get">
    <input type="text" name="q" placeholder="Search docs" />
    <input type="hidden" name="check_keywords" value="yes" />
    <input type="hidden" name="area" value="default" />
  </form>
</div>
        </div><div class="wy-menu wy-menu-vertical" data-spy="affix" role="navigation" aria-label="Navigation menu">
              <p class="caption" role="heading"><span class="caption-text">Contents:</span></p>
<ul class="current">
<li class="toctree-l1 current"><a class="reference internal" href="modules.html">bornraytrace</a><ul class="current">
<li class="toctree-l2 current"><a class="current reference internal" href="#">bornraytrace package</a><ul>
<li class="toctree-l3"><a class="reference internal" href="#submodules">Submodules</a></li>
<li class="toctree-l3"><a class="reference internal" href="#module-bornraytrace.intrinsic_alignments">bornraytrace.intrinsic_alignments module</a></li>
<li class="toctree-l3"><a class="reference internal" href="#module-bornraytrace.lensing">bornraytrace.lensing module</a></li>
<li class="toctree-l3"><a class="reference internal" href="#module-bornraytrace">Module contents</a></li>
</ul>
</li>
</ul>
</li>
</ul>

        </div>
      </div>
    </nav>

    <section data-toggle="wy-nav-shift" class="wy-nav-content-wrap"><nav class="wy-nav-top" aria-label="Mobile navigation menu" >
          <i data-toggle="wy-nav-top" class="fa fa-bars"></i>
          <a href="index.html">BornRaytrace</a>
      </nav>

      <div class="wy-nav-content">
        <div class="rst-content">
          <div role="navigation" aria-label="Page navigation">
  <ul class="wy-breadcrumbs">
      <li><a href="index.html" class="icon icon-home"></a> &raquo;</li>
          <li><a href="modules.html">bornraytrace</a> &raquo;</li>
      <li>bornraytrace package</li>
      <li class="wy-breadcrumbs-aside">
            <a href="_sources/bornraytrace.rst.txt" rel="nofollow"> View page source</a>
      </li>
  </ul>
  <hr/>
</div>
          <div role="main" class="document" itemscope="itemscope" itemtype="http://schema.org/Article">
           <div itemprop="articleBody">
             
  <section id="bornraytrace-package">
<h1>bornraytrace package<a class="headerlink" href="#bornraytrace-package" title="Permalink to this headline"></a></h1>
<section id="submodules">
<h2>Submodules<a class="headerlink" href="#submodules" title="Permalink to this headline"></a></h2>
</section>
<section id="module-bornraytrace.intrinsic_alignments">
<span id="bornraytrace-intrinsic-alignments-module"></span><h2>bornraytrace.intrinsic_alignments module<a class="headerlink" href="#module-bornraytrace.intrinsic_alignments" title="Permalink to this headline"></a></h2>
<dl class="py function">
<dt class="sig sig-object py" id="bornraytrace.intrinsic_alignments.D_1">
<span class="sig-prename descclassname"><span class="pre">bornraytrace.intrinsic_alignments.</span></span><span class="sig-name descname"><span class="pre">D_1</span></span><span class="sig-paren">(</span><em class="sig-param"><span class="n"><span class="pre">z</span></span></em>, <em class="sig-param"><span class="n"><span class="pre">om0</span></span></em><span class="sig-paren">)</span><a class="reference internal" href="_modules/bornraytrace/intrinsic_alignments.html#D_1"><span class="viewcode-link"><span class="pre">[source]</span></span></a><a class="headerlink" href="#bornraytrace.intrinsic_alignments.D_1" title="Permalink to this definition"></a></dt>
<dd><p>Normalised linear growth factor (D_plus)</p>
<dl class="field-list simple">
<dt class="field-odd">Parameters</dt>
<dd class="field-odd"><ul class="simple">
<li><p><strong>z</strong> – single redshift value or array values</p></li>
<li><p><strong>om0</strong> – matter density</p></li>
</ul>
</dd>
<dt class="field-even">Returns</dt>
<dd class="field-even"><p>normalised linear growth factor</p>
</dd>
</dl>
</dd></dl>

<dl class="py function">
<dt class="sig sig-object py" id="bornraytrace.intrinsic_alignments.D_single">
<span class="sig-prename descclassname"><span class="pre">bornraytrace.intrinsic_alignments.</span></span><span class="sig-name descname"><span class="pre">D_single</span></span><span class="sig-paren">(</span><em class="sig-param"><span class="n"><span class="pre">z</span></span></em>, <em class="sig-param"><span class="n"><span class="pre">om0</span></span></em><span class="sig-paren">)</span><a class="reference internal" href="_modules/bornraytrace/intrinsic_alignments.html#D_single"><span class="viewcode-link"><span class="pre">[source]</span></span></a><a class="headerlink" href="#bornraytrace.intrinsic_alignments.D_single" title="Permalink to this definition"></a></dt>
<dd><p>Provides the normalised linear growth factor</p>
<dl class="field-list simple">
<dt class="field-odd">Parameters</dt>
<dd class="field-odd"><ul class="simple">
<li><p><strong>z</strong> – single redshift value</p></li>
<li><p><strong>om0</strong> – matter density</p></li>
</ul>
</dd>
<dt class="field-even">Returns</dt>
<dd class="field-even"><p>normalised linear growth factor</p>
</dd>
</dl>
</dd></dl>

<dl class="py function">
<dt class="sig sig-object py" id="bornraytrace.intrinsic_alignments.E_sq">
<span class="sig-prename descclassname"><span class="pre">bornraytrace.intrinsic_alignments.</span></span><span class="sig-name descname"><span class="pre">E_sq</span></span><span class="sig-paren">(</span><em class="sig-param"><span class="n"><span class="pre">z</span></span></em>, <em class="sig-param"><span class="n"><span class="pre">om0</span></span></em><span class="sig-paren">)</span><a class="reference internal" href="_modules/bornraytrace/intrinsic_alignments.html#E_sq"><span class="viewcode-link"><span class="pre">[source]</span></span></a><a class="headerlink" href="#bornraytrace.intrinsic_alignments.E_sq" title="Permalink to this definition"></a></dt>
<dd><p>A function giving Hubble’s law for flat cosmology</p>
<dl class="field-list simple">
<dt class="field-odd">Parameters</dt>
<dd class="field-odd"><ul class="simple">
<li><p><strong>z</strong> – redshift value</p></li>
<li><p><strong>om0</strong> – matter density</p></li>
</ul>
</dd>
<dt class="field-even">Returns</dt>
<dd class="field-even"><p>A value for the Hubble parameter</p>
</dd>
</dl>
</dd></dl>

<dl class="py function">
<dt class="sig sig-object py" id="bornraytrace.intrinsic_alignments.F_nla">
<span class="sig-prename descclassname"><span class="pre">bornraytrace.intrinsic_alignments.</span></span><span class="sig-name descname"><span class="pre">F_nla</span></span><span class="sig-paren">(</span><em class="sig-param"><span class="n"><span class="pre">z</span></span></em>, <em class="sig-param"><span class="n"><span class="pre">om0</span></span></em>, <em class="sig-param"><span class="n"><span class="pre">A_ia</span></span></em>, <em class="sig-param"><span class="n"><span class="pre">rho_c1</span></span></em>, <em class="sig-param"><span class="n"><span class="pre">eta</span></span><span class="o"><span class="pre">=</span></span><span class="default_value"><span class="pre">0.0</span></span></em>, <em class="sig-param"><span class="n"><span class="pre">z0</span></span><span class="o"><span class="pre">=</span></span><span class="default_value"><span class="pre">0.0</span></span></em>, <em class="sig-param"><span class="n"><span class="pre">lbar</span></span><span class="o"><span class="pre">=</span></span><span class="default_value"><span class="pre">0.0</span></span></em>, <em class="sig-param"><span class="n"><span class="pre">l0</span></span><span class="o"><span class="pre">=</span></span><span class="default_value"><span class="pre">1e-09</span></span></em>, <em class="sig-param"><span class="n"><span class="pre">beta</span></span><span class="o"><span class="pre">=</span></span><span class="default_value"><span class="pre">0.0</span></span></em><span class="sig-paren">)</span><a class="reference internal" href="_modules/bornraytrace/intrinsic_alignments.html#F_nla"><span class="viewcode-link"><span class="pre">[source]</span></span></a><a class="headerlink" href="#bornraytrace.intrinsic_alignments.F_nla" title="Permalink to this definition"></a></dt>
<dd><p>NLA intrinsic alignment amplitude</p>
<dl class="field-list simple">
<dt class="field-odd">Parameters</dt>
<dd class="field-odd"><ul class="simple">
<li><p><strong>z</strong> – redshift value</p></li>
<li><p><strong>om0</strong> – matter density</p></li>
<li><p><strong>A_ia</strong> – amplitude parameter</p></li>
<li><p><strong>rho_c1</strong> – rho_crit x C1 (C1 approx 1.508e+27 cm3 / g)</p></li>
<li><p><strong>eta</strong> – redshift dependence</p></li>
<li><p><strong>z0</strong> – arbitrary redshift pivot parameter</p></li>
<li><p><strong>lbar</strong> – average luminosity of source galaxy population</p></li>
<li><p><strong>l0</strong> – arbitrary luminosity pivot parameter</p></li>
<li><p><strong>beta</strong> – luminosity dependence</p></li>
</ul>
</dd>
<dt class="field-even">Returns</dt>
<dd class="field-even"><p>NLA F(z) amplitude</p>
</dd>
</dl>
</dd></dl>

<dl class="py function">
<dt class="sig sig-object py" id="bornraytrace.intrinsic_alignments.f_integrand">
<span class="sig-prename descclassname"><span class="pre">bornraytrace.intrinsic_alignments.</span></span><span class="sig-name descname"><span class="pre">f_integrand</span></span><span class="sig-paren">(</span><em class="sig-param"><span class="n"><span class="pre">z</span></span></em>, <em class="sig-param"><span class="n"><span class="pre">om0</span></span></em><span class="sig-paren">)</span><a class="reference internal" href="_modules/bornraytrace/intrinsic_alignments.html#f_integrand"><span class="viewcode-link"><span class="pre">[source]</span></span></a><a class="headerlink" href="#bornraytrace.intrinsic_alignments.f_integrand" title="Permalink to this definition"></a></dt>
<dd><p>A function for the redshift integrand in the intrinsic alignment calculation</p>
<dl class="field-list simple">
<dt class="field-odd">Parameters</dt>
<dd class="field-odd"><ul class="simple">
<li><p><strong>z</strong> – redshift value</p></li>
<li><p><strong>om0</strong> – matter density</p></li>
</ul>
</dd>
<dt class="field-even">Returns</dt>
<dd class="field-even"><p>redshift integrand</p>
</dd>
</dl>
</dd></dl>

</section>
<section id="module-bornraytrace.lensing">
<span id="bornraytrace-lensing-module"></span><h2>bornraytrace.lensing module<a class="headerlink" href="#module-bornraytrace.lensing" title="Permalink to this headline"></a></h2>
<dl class="py function">
<dt class="sig sig-object py" id="bornraytrace.lensing.W_kernel">
<span class="sig-prename descclassname"><span class="pre">bornraytrace.lensing.</span></span><span class="sig-name descname"><span class="pre">W_kernel</span></span><span class="sig-paren">(</span><em class="sig-param"><span class="n"><span class="pre">r_array</span></span></em>, <em class="sig-param"><span class="n"><span class="pre">z_array</span></span></em>, <em class="sig-param"><span class="n"><span class="pre">nz</span></span></em>, <em class="sig-param"><span class="n"><span class="pre">simpsons</span></span><span class="o"><span class="pre">=</span></span><span class="default_value"><span class="pre">False</span></span></em><span class="sig-paren">)</span><a class="reference internal" href="_modules/bornraytrace/lensing.html#W_kernel"><span class="viewcode-link"><span class="pre">[source]</span></span></a><a class="headerlink" href="#bornraytrace.lensing.W_kernel" title="Permalink to this definition"></a></dt>
<dd><p>lensing kernel W s.t.  kappa = prefactor * integral  W(r) * overdensity(r)  dr</p>
<dl class="field-list simple">
<dt class="field-odd">Parameters</dt>
<dd class="field-odd"><ul class="simple">
<li><p><strong>r_array</strong> – comoving distances array</p></li>
<li><p><strong>z_array</strong> – redshift array matching r_array (cosmology dependent)</p></li>
<li><p><strong>nz</strong> – source redshift distribution</p></li>
<li><p><strong>simpsons</strong> – boolean to use simpsons integratio</p></li>
</ul>
</dd>
<dt class="field-even">Returns</dt>
<dd class="field-even"><p>W = r * q /r</p>
</dd>
</dl>
</dd></dl>

<dl class="py function">
<dt class="sig sig-object py" id="bornraytrace.lensing.get_neighbour_array">
<span class="sig-prename descclassname"><span class="pre">bornraytrace.lensing.</span></span><span class="sig-name descname"><span class="pre">get_neighbour_array</span></span><span class="sig-paren">(</span><em class="sig-param"><span class="n"><span class="pre">nside</span></span></em><span class="sig-paren">)</span><a class="reference internal" href="_modules/bornraytrace/lensing.html#get_neighbour_array"><span class="viewcode-link"><span class="pre">[source]</span></span></a><a class="headerlink" href="#bornraytrace.lensing.get_neighbour_array" title="Permalink to this definition"></a></dt>
<dd><p>array of indices labelling the 8 neighbouring pixels for each pixel</p>
<dl class="field-list simple">
<dt class="field-odd">Parameters</dt>
<dd class="field-odd"><p><strong>nside</strong> – nside of map</p>
</dd>
<dt class="field-even">Returns</dt>
<dd class="field-even"><p>neighbour indices array</p>
</dd>
</dl>
</dd></dl>

<dl class="py function">
<dt class="sig sig-object py" id="bornraytrace.lensing.kappa2shear">
<span class="sig-prename descclassname"><span class="pre">bornraytrace.lensing.</span></span><span class="sig-name descname"><span class="pre">kappa2shear</span></span><span class="sig-paren">(</span><em class="sig-param"><span class="n"><span class="pre">kappa_map</span></span></em>, <em class="sig-param"><span class="n"><span class="pre">lmax</span></span><span class="o"><span class="pre">=</span></span><span class="default_value"><span class="pre">None</span></span></em><span class="sig-paren">)</span><a class="reference internal" href="_modules/bornraytrace/lensing.html#kappa2shear"><span class="viewcode-link"><span class="pre">[source]</span></span></a><a class="headerlink" href="#bornraytrace.lensing.kappa2shear" title="Permalink to this definition"></a></dt>
<dd><p>Performs inverse Kaiser-Squires on the sphere with healpy spherical harmonics</p>
<dl class="field-list simple">
<dt class="field-odd">Parameters</dt>
<dd class="field-odd"><ul class="simple">
<li><p><strong>kappa_map</strong> – healpix format complex convergence (kappa) map</p></li>
<li><p><strong>lmax</strong> – maximum multipole</p></li>
</ul>
</dd>
<dt class="field-even">Returns</dt>
<dd class="field-even"><p>complex shear map (gamma1 + 1j * gamma2)</p>
</dd>
</dl>
</dd></dl>

<dl class="py function">
<dt class="sig sig-object py" id="bornraytrace.lensing.kappa_prefactor">
<span class="sig-prename descclassname"><span class="pre">bornraytrace.lensing.</span></span><span class="sig-name descname"><span class="pre">kappa_prefactor</span></span><span class="sig-paren">(</span><em class="sig-param"><span class="n"><span class="pre">H0</span></span></em>, <em class="sig-param"><span class="n"><span class="pre">om0</span></span></em>, <em class="sig-param"><span class="n"><span class="pre">length_unit</span></span><span class="o"><span class="pre">=</span></span><span class="default_value"><span class="pre">'Mpc'</span></span></em><span class="sig-paren">)</span><a class="reference internal" href="_modules/bornraytrace/lensing.html#kappa_prefactor"><span class="viewcode-link"><span class="pre">[source]</span></span></a><a class="headerlink" href="#bornraytrace.lensing.kappa_prefactor" title="Permalink to this definition"></a></dt>
<dd><p>Gives prefactor (3 H_0^2 Om0)/2</p>
<dl class="field-list simple">
<dt class="field-odd">Parameters</dt>
<dd class="field-odd"><ul class="simple">
<li><p><strong>H0</strong> – Hubble parameter with astropy units</p></li>
<li><p><strong>om0</strong> – Omega matter</p></li>
<li><p><strong>length_unit</strong> – for H0 (default Mpc)</p></li>
</ul>
</dd>
<dt class="field-even">Returns</dt>
<dd class="field-even"><p>prefactor for lensing</p>
</dd>
</dl>
</dd></dl>

<dl class="py function">
<dt class="sig sig-object py" id="bornraytrace.lensing.peak_find">
<span class="sig-prename descclassname"><span class="pre">bornraytrace.lensing.</span></span><span class="sig-name descname"><span class="pre">peak_find</span></span><span class="sig-paren">(</span><em class="sig-param"><span class="n"><span class="pre">map_input</span></span></em>, <em class="sig-param"><span class="n"><span class="pre">nside</span></span></em>, <em class="sig-param"><span class="n"><span class="pre">neighbour_array</span></span><span class="o"><span class="pre">=</span></span><span class="default_value"><span class="pre">None</span></span></em><span class="sig-paren">)</span><a class="reference internal" href="_modules/bornraytrace/lensing.html#peak_find"><span class="viewcode-link"><span class="pre">[source]</span></span></a><a class="headerlink" href="#bornraytrace.lensing.peak_find" title="Permalink to this definition"></a></dt>
<dd><p>Find peaks (local maxima) for a given input map</p>
<dl class="field-list simple">
<dt class="field-odd">Parameters</dt>
<dd class="field-odd"><ul class="simple">
<li><p><strong>map_input</strong> – input map</p></li>
<li><p><strong>nside</strong> – nside of map</p></li>
<li><p><strong>neighbour_array</strong> – optional array of indices labelling the 8 neighbouring pixels for each pixel</p></li>
</ul>
</dd>
<dt class="field-even">Returns</dt>
<dd class="field-even"><p>list of pixel indices for the peaks</p>
</dd>
</dl>
</dd></dl>

<dl class="py function">
<dt class="sig sig-object py" id="bornraytrace.lensing.raytrace">
<span class="sig-prename descclassname"><span class="pre">bornraytrace.lensing.</span></span><span class="sig-name descname"><span class="pre">raytrace</span></span><span class="sig-paren">(</span><em class="sig-param"><span class="n"><span class="pre">H0</span></span></em>, <em class="sig-param"><span class="n"><span class="pre">om0</span></span></em>, <em class="sig-param"><span class="n"><span class="pre">overdensity_array</span></span></em>, <em class="sig-param"><span class="n"><span class="pre">a_centre</span></span></em>, <em class="sig-param"><span class="n"><span class="pre">comoving_edges</span></span></em>, <em class="sig-param"><span class="n"><span class="pre">mask</span></span><span class="o"><span class="pre">=</span></span><span class="default_value"><span class="pre">None</span></span></em>, <em class="sig-param"><span class="n"><span class="pre">Hubble_length_unit</span></span><span class="o"><span class="pre">=</span></span><span class="default_value"><span class="pre">'Mpc'</span></span></em><span class="sig-paren">)</span><a class="reference internal" href="_modules/bornraytrace/lensing.html#raytrace"><span class="viewcode-link"><span class="pre">[source]</span></span></a><a class="headerlink" href="#bornraytrace.lensing.raytrace" title="Permalink to this definition"></a></dt>
<dd><p>Evaluate weak lensing convergence map using Born approximation</p>
<dl class="field-list simple">
<dt class="field-odd">Parameters</dt>
<dd class="field-odd"><ul class="simple">
<li><p><strong>H0</strong> – Hubble parameter with astropy units</p></li>
<li><p><strong>om0</strong> – Omega matter</p></li>
<li><p><strong>overdensity_array</strong> – an 2D array of overdensity healpix maps in radial shells</p></li>
<li><p><strong>a_centre</strong> – scale factor at comoving centre of shells</p></li>
<li><p><strong>comoving_edges</strong> – comoving distance to edges of shells</p></li>
<li><p><strong>mask</strong> – healpix map where 1 is observed and 0 is mask</p></li>
<li><p><strong>length_unit</strong> – for H0 (default Mpc)</p></li>
</ul>
</dd>
<dt class="field-even">Returns</dt>
<dd class="field-even"><p>convergence kappa map</p>
</dd>
</dl>
</dd></dl>

<dl class="py function">
<dt class="sig sig-object py" id="bornraytrace.lensing.raytrace_integration">
<span class="sig-prename descclassname"><span class="pre">bornraytrace.lensing.</span></span><span class="sig-name descname"><span class="pre">raytrace_integration</span></span><span class="sig-paren">(</span><em class="sig-param"><span class="n"><span class="pre">kappa_prefactor</span></span></em>, <em class="sig-param"><span class="n"><span class="pre">overdensity_array</span></span></em>, <em class="sig-param"><span class="n"><span class="pre">a_centre</span></span></em>, <em class="sig-param"><span class="n"><span class="pre">comoving_edges</span></span></em>, <em class="sig-param"><span class="n"><span class="pre">mask</span></span><span class="o"><span class="pre">=</span></span><span class="default_value"><span class="pre">None</span></span></em><span class="sig-paren">)</span><a class="reference internal" href="_modules/bornraytrace/lensing.html#raytrace_integration"><span class="viewcode-link"><span class="pre">[source]</span></span></a><a class="headerlink" href="#bornraytrace.lensing.raytrace_integration" title="Permalink to this definition"></a></dt>
<dd><p>This function evaluates the Born weak lensing integral</p>
<dl class="field-list simple">
<dt class="field-odd">Parameters</dt>
<dd class="field-odd"><ul class="simple">
<li><p><strong>kappa_prefactor</strong> – defined as the output of the function kappa_prefactor</p></li>
<li><p><strong>overdensity_array</strong> – an 2D array of overdensity healpix maps in radial shells</p></li>
<li><p><strong>a_centre</strong> – scale factor at comoving centre of shells</p></li>
<li><p><strong>comoving_edges</strong> – comoving distance to edges of shells</p></li>
<li><p><strong>mask</strong> – healpix map where 1 is observed and 0 is mask</p></li>
</ul>
</dd>
<dt class="field-even">Returns</dt>
<dd class="field-even"><p>convergence kappa map</p>
</dd>
</dl>
</dd></dl>

<dl class="py function">
<dt class="sig sig-object py" id="bornraytrace.lensing.recentre_nz">
<span class="sig-prename descclassname"><span class="pre">bornraytrace.lensing.</span></span><span class="sig-name descname"><span class="pre">recentre_nz</span></span><span class="sig-paren">(</span><em class="sig-param"><span class="n"><span class="pre">z_sim_edges</span></span></em>, <em class="sig-param"><span class="n"><span class="pre">z_samp_centre</span></span></em>, <em class="sig-param"><span class="n"><span class="pre">nz_input</span></span></em><span class="sig-paren">)</span><a class="reference internal" href="_modules/bornraytrace/lensing.html#recentre_nz"><span class="viewcode-link"><span class="pre">[source]</span></span></a><a class="headerlink" href="#bornraytrace.lensing.recentre_nz" title="Permalink to this definition"></a></dt>
<dd><p>Takes input n(z) sampled at z_samp_centre
and evaluates interpolated n(z) at new z values
to match a simulation at z_sim_edges</p>
<dl class="field-list simple">
<dt class="field-odd">Parameters</dt>
<dd class="field-odd"><ul class="simple">
<li><p><strong>z_sim_edges</strong> – new z values for n(z)</p></li>
<li><p><strong>z_samp_centre</strong> – original z values for n(z)</p></li>
<li><p><strong>nz_input</strong> – original n(z)</p></li>
</ul>
</dd>
<dt class="field-even">Returns</dt>
<dd class="field-even"><p>new n(z)</p>
</dd>
</dl>
</dd></dl>

<dl class="py function">
<dt class="sig sig-object py" id="bornraytrace.lensing.rotate_mask_approx">
<span class="sig-prename descclassname"><span class="pre">bornraytrace.lensing.</span></span><span class="sig-name descname"><span class="pre">rotate_mask_approx</span></span><span class="sig-paren">(</span><em class="sig-param"><span class="n"><span class="pre">mask</span></span></em>, <em class="sig-param"><span class="n"><span class="pre">rot_angles</span></span></em>, <em class="sig-param"><span class="n"><span class="pre">flip</span></span><span class="o"><span class="pre">=</span></span><span class="default_value"><span class="pre">False</span></span></em><span class="sig-paren">)</span><a class="reference internal" href="_modules/bornraytrace/lensing.html#rotate_mask_approx"><span class="viewcode-link"><span class="pre">[source]</span></span></a><a class="headerlink" href="#bornraytrace.lensing.rotate_mask_approx" title="Permalink to this definition"></a></dt>
<dd><p>rotate healpix mask on sphere</p>
<dl class="field-list simple">
<dt class="field-odd">Parameters</dt>
<dd class="field-odd"><ul class="simple">
<li><p><strong>mask</strong> – healpix map of ones and zeros</p></li>
<li><p><strong>rot_angles</strong> – rotation on the sphere (e.g. [ 45.91405291 ,150.72092269 , 46.34505909])</p></li>
<li><p><strong>flip</strong> – boolean, mirror the mask</p></li>
</ul>
</dd>
<dt class="field-even">Returns</dt>
<dd class="field-even"><p>rotated map</p>
</dd>
</dl>
</dd></dl>

<dl class="py function">
<dt class="sig sig-object py" id="bornraytrace.lensing.shear2kappa">
<span class="sig-prename descclassname"><span class="pre">bornraytrace.lensing.</span></span><span class="sig-name descname"><span class="pre">shear2kappa</span></span><span class="sig-paren">(</span><em class="sig-param"><span class="n"><span class="pre">shear_map</span></span></em>, <em class="sig-param"><span class="n"><span class="pre">lmax</span></span><span class="o"><span class="pre">=</span></span><span class="default_value"><span class="pre">None</span></span></em><span class="sig-paren">)</span><a class="reference internal" href="_modules/bornraytrace/lensing.html#shear2kappa"><span class="viewcode-link"><span class="pre">[source]</span></span></a><a class="headerlink" href="#bornraytrace.lensing.shear2kappa" title="Permalink to this definition"></a></dt>
<dd><p>Performs Kaiser-Squires on the sphere with healpy spherical harmonics</p>
<dl class="field-list simple">
<dt class="field-odd">Parameters</dt>
<dd class="field-odd"><ul class="simple">
<li><p><strong>shear_map</strong> – healpix format complex shear map</p></li>
<li><p><strong>lmax</strong> – maximum ell multipole</p></li>
</ul>
</dd>
<dt class="field-even">Returns</dt>
<dd class="field-even"><p>kappa map</p>
</dd>
</dl>
</dd></dl>

</section>
<section id="module-bornraytrace">
<span id="module-contents"></span><h2>Module contents<a class="headerlink" href="#module-bornraytrace" title="Permalink to this headline"></a></h2>
</section>
</section>


           </div>
          </div>
          <footer><div class="rst-footer-buttons" role="navigation" aria-label="Footer">
        <a href="modules.html" class="btn btn-neutral float-left" title="bornraytrace" accesskey="p" rel="prev"><span class="fa fa-arrow-circle-left" aria-hidden="true"></span> Previous</a>
    </div>

  <hr/>

  <div role="contentinfo">
    <p>&#169; Copyright 2021, Niall Jeffrey.</p>
  </div>

  Built with <a href="https://www.sphinx-doc.org/">Sphinx</a> using a
    <a href="https://github.com/readthedocs/sphinx_rtd_theme">theme</a>
    provided by <a href="https://readthedocs.org">Read the Docs</a>.
   

</footer>
        </div>
      </div>
    </section>
  </div>
  <script>
      jQuery(function () {
          SphinxRtdTheme.Navigation.enable(true);
      });
  </script> 

</body>
</html>

# eraser_afshar
Simulation of the quantuim eraser plus Afshar's experiment (and H and V polarizers in front of the detectors)

The quantum eraser experiment used one beam of entangled photons to control interference in another beam of
entangled partner photons.  The original experiment used coincidence detection to filter one of the two
interference patterns - they would overlap, washing out the interference.  However, by adding Afshar's
experiment to the quantum eraser, coincidence is no longer needed.

This is a simulation where light of a certain polarization is sent through two slits with quarter wave plates
in front of them.  There are beam block, a lens, and detectors.  This has three files: the main Python sim,
an object to handle calculations for interference and ray tracing, and a JSON file to configure params of the
optical setup, what light is sent in, etc.

To run, use python to execute sim2.py and give it a name from the config file (like "first").  Then install
the packages that it needs and try again :)

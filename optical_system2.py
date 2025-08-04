import sys
import numpy as np
import math


# class for optical setup

class OpticalSystem:
    """
    Holds parameters for system and provides functions for calculating interference, etc.
    """

    # polarizations

    pol_H = np.array([1, 0], dtype=complex)
    pol_V = np.array([0, 1], dtype=complex)

    pol_D = np.array([1, 1], dtype=complex) / np.sqrt(2)
    pol_A = np.array([1, -1], dtype=complex) / np.sqrt(2)


    # translate betwwen degrees and radians

    dtor = (2*np.pi / 360.0)
    rtod = (360.0 / 2*np.pi)


    def __init__(self, qwp_upper=None, qwp_lower=None, wavelength=0.1,
                 x_lens=50.0, r_lens=100, n_lens=1.3,
                 x_pos=3.0, y_delta=2.0,
                 field_width=100, field_height=100,
                 scale_factor=1.0):

        """
        Initialize optical system
        inputs:
            upper quarter wave plate angle (rotation degrees or None)
            lower quarter wave plate angle (rotation degrees or None)
            x position of lens (mm)
            x position of light sources (mm)
            y delta between sources (mm)
            field width (mm)
            field height (mm)
            scale factor for angle intensity falloff (float)
        """

        # create QWP objects if angles are specified

        if qwp_upper is not None and qwp_upper >= 0:
            self.qwp_upper = self.qwp(qwp_upper)
        else:
            self.qwp_upper = None

        if qwp_lower is not None and qwp_lower >= 0:
            self.qwp_lower = self.qwp(qwp_lower)
        else:
            self.qwp_lower = None


        # save wavelength of light

        self.wavelength = wavelength


        # save positions of lens and sources

        self.x_lens = x_lens
        self.r_lens = r_lens
        self.n_lens = n_lens

        self.x_sources = x_pos
        self.y_source_pos = y_delta / 2.0


        # save field size

        self.x_size = field_width
        self.y_size = field_height


        # intensity falloff scale factor

        self.scale_factor = scale_factor


    def qwp(self, theta_deg):
        # create rotation matrix

        theta = theta_deg * self.dtor

        rot_p = np.array([
            [np.cos(theta), -np.sin(theta)],
            [np.sin(theta), np.cos(theta)]
        ], dtype=complex)

        rot_m = np.array([
            [np.cos(-theta), -np.sin(-theta)],
            [np.sin(-theta), np.cos(-theta)]
        ], dtype=complex)


        # create QWP matrix

        qwp_basic = np.array([
            [1, 0],
            [0, 1j]
        ], dtype=complex)


        # return matrix

        return rot_m @ qwp_basic @ rot_p


    def angle_falloff(self, source_pos, point_pos):
        """
        Compute angle-based amplitude falloff (e.g., cos^2 of emission angle).
        Assumes emission is strongest along source's forward axis (default: +x).
        """

        # calculae vector from source to point

        vec = point_pos - source_pos             # from source to point

        if np.linalg.norm(vec) == 0:
            return 0.0  # Avoid divide by zero

        unit_vec = vec / np.linalg.norm(vec)


        # make horizontal unit vector, calculate angle and scale

        horizontal = np.array([1.0, 0.0])

        theta = np.arccos(np.clip(np.dot(unit_vec, horizontal), -1.0, 1.0))
        scaled_theta = self.scale_factor * theta

        if scaled_theta > np.pi/2:
            return 0.0  # greater than 90 degrees


        # return scale factor

        v = np.cos(scaled_theta)
        return v


    def distance_falloff(self, source_pos, point_pos):
        """
        Compute amplitude falloff due to distance (e.g., 1/r).
        You could use 1/r or 1/r^2 depending on physical assumptions.
        """

        r = np.linalg.norm(point_pos - source_pos)

        if r == 0:
            return 1.0  # source at point

        return 1.0 / r


    def interference_factor(self, pol, x, y):
        """
        Calculate interference factor at a point
        """

        # source positions

        pos_upper = np.array([self.x_sources, -0.5 * self.y_source_pos])
        pos_lower = np.array([self.x_sources, 0.5 * self.y_source_pos])

        point = np.array([x, y])


        # distances from sources to point

        r_upper = np.linalg.norm(point - pos_upper)
        r_lower = np.linalg.norm(point - pos_lower)


        # phase from each source

        phase_upper = np.exp(1j * 2 * np.pi * r_upper / self.wavelength)
        phase_lower = np.exp(1j * 2 * np.pi * r_lower / self.wavelength)


        # Amplitudes due to angle and distance falloff

        A_upper = self.angle_falloff(pos_upper, point) * self.distance_falloff(pos_upper, point)
        A_lower = self.angle_falloff(pos_lower, point) * self.distance_falloff(pos_lower, point)


        # create E-field ectors and apply QWPs

        E_upper = pol * A_upper
        E_lower = pol * A_lower

        if self.qwp_upper is not None:
            E_upper = self.qwp_upper @ E_upper

        if self.qwp_lower is not None:
            E_lower = self.qwp_lower @ E_lower


        # apply phase

        E_upper *= phase_upper
        E_lower *= phase_lower


        # computer intensity

        E_total = E_upper + E_lower
        I_total = np.vdot(E_total, E_total).real


        # compute max intensity

        I_upper = np.vdot(E_upper, E_upper).real
        I_lower = np.vdot(E_lower, E_lower).real

        I_max = I_upper + I_lower + 2*np.sqrt(I_upper * I_lower)


        # return interference factor

        if I_max == 0:
            return 0.0

        else:
            return I_total / I_max


    def lens_angle(self, offset, angle_in):
        # lens

        x_lens = self.x_lens

        r_lens = self.r_lens
        n_lens = self.n_lens


        # Surface angle (normal to curved surface)

        if offset/r_lens > 1.0 or offset/r_lens < -1.0:
            return angle_in

        surface_angle = math.asin(offset/r_lens)

        # Incident angle relative to surface normal
        incident_relative = angle_in - surface_angle

        # Snell's law at first surface: air -> glass
        refracted_relative = math.asin((1/n_lens) * math.sin(incident_relative))

        # For thin symmetric lens, second surface gives opposite deflection
        # Exit angle relative to optical axis
        exit_angle = angle_in + 2 * (incident_relative - refracted_relative)

        return exit_angle


    def ray_trace(self, polarization, x_initial, y_initial, ray_angle):
        # build up segments

        ray_points = []
        ray_points.append([x_initial, y_initial])

        ray_continue = True


        # check for stupid angles

        if ray_angle == 0.0:
            x_int = self.x_lens
            y_int = y_initial

            ray_points.append([x_int, y_int])

        elif ray_angle == -90.0:
            x_int = x_initial
            y_int = -field_height/2

            ray_points.append([x_int, y_int])
            ray_continue = False

        elif ray_angle == 90.0:
            x_int = x_initial
            y_int = field_height/2

            ray_points.append([x_int, y_int])
            ray_continue = False

        else:
            # calc slope, find intersection with lens

            ray_slope = math.tan(ray_angle * self.dtor)
            y_int = y_initial + self.x_lens * ray_slope

            if y_int < -self.y_size/2:
                y_int = -self.y_size/2
                x_int = y_int / ray_slope

                ray_points.append([x_int, y_int])
                ray_continue = False

            elif y_int > self.y_size/2:
                y_int = self.y_size/2
                x_int = y_int / ray_slope

                ray_points.append([x_int, y_int])
                ray_continue = False

            else:
                x_int = self.x_lens

                ray_points.append([x_int, y_int])


        # calc intensity

        inten_interferencee = self.interference_factor(polarization, x_int, y_int)

        pos_source = np.array([x_initial, y_initial])
        point = np.array([x_int, y_int])

        inten_angle = self.angle_falloff(pos_source, point)


        # if ray hit lens, then continue ray

        if ray_continue:
            # calculate new angle

            ray_angle2 = self.lens_angle(y_int, ray_angle * self.dtor)


            # and again

            if ray_angle2 == 0.0:
                x_int2 = self.x_size
                y_int2 = y_int

                ray_points.append((x_int2, y_int2))

            else:
                # get slope of ray, intersect with wall

                ray_slope2 = math.tan(ray_angle2)
                y_int2 = y_int + (self.x_size - self.x_lens) * ray_slope2

                if y_int2 < -self.y_size/2:
                    y_int2 = -self.y_size/2
                    x_int2 = y_int2 / ray_slope2

                    ray_points.append([x_int2, y_int2])

                elif y_int2 > self.y_size/2:
                    y_int2 = self.y_size/2
                    x_int2 = y_int2 / ray_slope2

                    ray_points.append([x_int2, y_int2])

                else:
                    x_int2 = self.x_size

                    ray_points.append([x_int2, y_int2])


        # return array of points

        return [ray_points, inten_interferencee*inten_angle]


    def trace_source(self, polarization, upper_lower, angle):
        if upper_lower == "upper":
            return self.ray_trace(polarization, self.x_sources, self.y_source_pos, angle)

        elif upper_lower == "lower":
            return self.ray_trace(polarization, self.x_sources, -self.y_source_pos, angle)

        else:
            print("What the cheese sandwich is this: " + str(upper_lower))
            sys.exit(0)

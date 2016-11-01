package G2;

use strict;
use Carp;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK $AUTOLOAD);

require Exporter;
require DynaLoader;
require AutoLoader;
use strict;

@ISA = qw(Exporter DynaLoader);
# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.
@EXPORT = qw(
	G2LD
	G2_H
	G2_VERSION
);
$VERSION = '0.01';

sub AUTOLOAD {
    # This AUTOLOAD is used to 'autoload' constants from the constant()
    # XS function.  If a constant is not found then control is passed
    # to the AUTOLOAD in AutoLoader.

    my $constname;
    ($constname = $AUTOLOAD) =~ s/.*:://;
    my $val = constant($constname, @_ ? $_[0] : 0);
    if ($! != 0) {
	if ($! =~ /Invalid/) {
	    $AutoLoader::AUTOLOAD = $AUTOLOAD;
	    goto &AutoLoader::AUTOLOAD;
	}
	else {
		croak "Your vendor has not defined G2 macro $constname";
	}
    }
    eval "sub $AUTOLOAD { $val }";
    goto &$AUTOLOAD;
}

bootstrap G2 $VERSION;

# Preloaded methods go here.

# Autoload methods go after =cut, and are processed by the autosplit program.

1;
__END__
# Below is the stub of documentation for your module. You better edit it!

=head1 NAME

G2 - A simple graphics library ported to Perl. 

=head1 SYNOPSIS

  use G2;

  $dev1 = newX11 G2::Device(775, 575);
  $dev2 = newGD  G2::Device("test.png",600,200);

  $dev1->rectangle(20,20,150,150);
  $dev1->circle(100,150,60);

  $dev2->circle(100,150,60);
  $dev2->string(100,50,"A circle in a PNG file");

=head1 DESCRIPTION

g2 is a simple to use graphics library for 2D graphical applications.
This library provides a comprehensive set of functions for
simultaneous generation of graphical output on different types of devices.
Presently, following devices are currently supported by g2: X11, PNG,
PostScript (xfig is in developement).
One major feature of the g2_library is the concept of virtual devices. An
arbitrary number of physical devices (such as PNG, or X11) can be grouped to
create a so-called virtual device. Commands sent to such a virtual devices
will automatically issued to all attached physical devices. This allows for
example simultaneous output to a PNG file and a Postscript file. A virtual
device in turn can be attached to another virtual device, allowing to
construct trees of devices.
Virtual devices can also be useful when using different user-coordinate
systems. E.g. one X11 window showing an overview of a graphical output, and
a second window showing  a zoom of a more detailed area of the graphic.
Drawing in both windows is performed by one single command to the virtual
device.
Please see g2 documentation (C interface) for up to date version.

=head1 Exported constants

  G2LD
  G2_H
  G2_VERSION


=head1 Exported functions

=head2 Creating new devices

=over 5

=item C<newX11>

C<G2::Device::newX11(width,height)> I<class method>

opens an X11 window with width and height of X11 window given in pixels.
returns : a new X11 device.

=item C<newPS>

C<G2::Device::newPS(filename,paper,orientation)> I<class method>

opens a new PostScript device.
file_name: name of PostScript file
paper: Paper size (e.g. g2_A4, g2_Letter). See
PostScript paper sizes for a full list of supported sizes.
orientation: paper orientation. Either g2_PS_land for
landscape or g2_PS_port  for  portrait
returns : a new PostScript device.

=item C<newGD>

C<G2::Device::newGD(filename,width,height,type)> I<class method>

open a new GD device
width,height: width and height of the image in pixels
filename: name of the output file.
type: file type, 0-jpeg, 1-png
returns : a new GD device

=item C<newVD>

C<G2::Device::newVD()> I<class method>

Create a new Virtual Device.  An arbitrary number of physical devices
(such as PNG, or X11) can be grouped to
create a so-called virtual device. Commands sent to such a virtual devices
will automatically issued to all attached physical devices. This allows for
example simultaneous output to a PNG file and a Postscript file. A virtual
device in turn can be attached to another virtual device, allowing to
construct trees of devices.
Virtual devices can also be useful when using different user-coordinate
systems. E.g. one X11 window showing an overview of a graphical output, and
a second window showing  a zoom of a more detailed area of the graphic.
Drawing in both windows is performed by one single command to the virtual
device.

=head2 Device Functions

=item C<>

C<G2::Device::attach(dev)> I<object method>

=item C<>

C<G2::Device::detach(dev)> I<object method>

=item C<>

C<G2::Device::close()> I<object method>

=item C<set_auto_flush>

C<G2::Device::set_auto_flush(on_off)> I<object method>

=item C<flush>

C<G2::Device::flush()> I<object method>

=item C<save>

C<G2::Device::save()> I<object method>

=item C<set_coordinate_system>

C<G2::Device::set_coordinate_system(x_origin, y_origin, x_mul, y_mul)> I<object method>

=item C<ld>

C<G2::Device::ld()> I<object method>

=item C<set_ld>

C<G2::Device::set_ld( dev)> I<object method>

=item C<ink>

C<G2::Device::ink( pd_dev, red, green, blue)> I<object method>

=item C<pen>

C<G2::Device::pen(color)> I<object method>

=item C<set_dash>

C<G2::Device::set_dash(N, *dashes)> I<object method>

=item C<set_font_size>

C<G2::Device::set_font_size(size)> I<object method>

=item C<set_line_width>

C<G2::Device::set_line_width(w)> I<object method>

=item C<clear_palette>

C<G2::Device::clear_palette( dev)> I<object method>

=item C<reset_palette>

C<G2::Device::reset_palette( dev)> I<object method>

=item C<allocate_basic_colors>

C<G2::Device::allocate_basic_colors( dev)> I<object method>

=item C<clear>

C<G2::Device::clear( dev)> I<object method>

=item C<set_background>

C<G2::Device::set_background(color)> I<object method>

=head2 Drawing Functions


=item C<move>

C<G2::Device::move(x, y)> I<object method>

=item C<move_r>

C<G2::Device::move_r(dx, dy)> I<object method>

=item C<plot>

C<G2::Device::plot(x, y)> I<object method>

=item C<plot_r>

C<G2::Device::plot_r(dx, dy)> I<object method>

=item C<line>

C<G2::Device::line(x1, y1, x2, y2)> I<object method>

=item C<line_r>

C<G2::Device::line_r(dx, dy)> I<object method>

=item C<line_to>

C<G2::Device::line_to(x, y)> I<object method>

=item C<poly_line>

C<G2::Device::poly_line(N_pt, *pos)> I<object method>

=item C<triangle>

C<G2::Device::triangle(x1, y1, x2, y2, x3, y3)> I<object method>

=item C<filled_triangle>

C<G2::Device::filled_triangle(x1, y1, x2, y2, x3, y3)> I<object method>

=item C<rectangle>


C<G2::Device::rectangle(x1, y1, x2, y2)> I<object method>

=item C<filled_rectangle>

C<G2::Device::filled_rectangle(x1, y1, x2, y2)> I<object method>

=item C<polygon>

C<G2::Device::polygon(N_pt, *pos)> I<object method>

=item C<filled_polygon>

C<G2::Device::filled_polygon(N_pt, *pos)> I<object method>

=item C<circle>

C<G2::Device::circle(x, y, r)> I<object method>

=item C<filled_circle>
	
C<G2::Device::filled_circle(x, y, r)> I<object method>

=item C<ellipse>

C<G2::Device::ellipse(x, y, r1, r2)> I<object method>
          Draw an ellipse on device dev
          x,y: center point
          r1,r2: x and y radius

=item C<filled_ellipse>

C<G2::Device::filled_ellipse(x, y, r1, r2)> I<object method>
          Draw a filled ellipse on device dev
          x,y: center point
          r1,r2: x and y radius

=item C<arc>

C<G2::Device::arc(x, y, r1, r2, a1, a2)> I<object method>

Draw an arc with center point at (x,y), x and y radius given by
r1,r2 and starting and ending angle in radians a1,a2

=item C<filled_arc>

C<G2::Device::filled_arc(x, y, r1, r2, a1, a2)> I<object method>
          Draw a filled arc on device dev
          x,y: center point
          r1,r2: x and y radius
          a1,a2: starting and ending angle in radians

=item C<string>

C<G2::Device::string(x, y, char *text)> I<object method>

=item C<>

C<G2::Device::set_QP(d, enum QPshape shape)> I<object method>

=item C<>

C<G2::Device::plot_QP(x, y)> I<object method>

=head1 AUTHORS

Horst Wagner (wagner/users-sourceforge.net) and Ljubomir Milanovic (ljubo/users-sourceforge-net)

=head1 COPYRIGHT

Copyright (C) 1998-2001  Ljubomir Milanovic & Horst Wagner
This file is part of the g2 library

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
    
=cut

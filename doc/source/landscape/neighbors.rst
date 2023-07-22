Neighborhood Relation and Move Sets for Secondary Structures
============================================================

Different functions to generate structural neighbors of a secondary
structure according to a particular Move Set.

This module contains methods to compute the neighbors of an RNA
secondary structure. Neighbors of a given structure are all structures
that differ in exactly one base pair. That means one can insert an
delete base pairs in the given structure. These insertions and
deletions of base pairs are usually called moves. A third move which
is considered in these methods is a shift move. A shifted base pair
has one stable position and one position that changes. These moves
are encoded as follows:

- insertion: :math:`(i, j)` where :math:`i,j > 0`
- deletion: :math:`(i, j)` where :math:`i,j` < 0
- shift: :math:`(i, j)` where either :math:`i > 0, j < 0` or :math:`i < 0, j > 0`

The negative position of a shift indicates the position that has changed.

An example::

    We are given a sequence and a structure.
    Sequence  AAGGAAACC
     Structure ..(.....)
     Indices   123456789

     The given base pair is (3,9) and the neighbors are the insertion (4, 8),
     the deletion (-3,-9), the shift (3,-8) and the shift (-4, 9).
     This leads to the neighbored structures:
     ...(....)
     .........
     ...(...).
     ....(...)

A simple method to construct all insertions is to iterate over the positions of a sequence twice. The first
iteration has the index i in [1, sequence length], the second iteration has the index j in [i+1, sequence length].
All pairs (i,j) with compatible letters and which are non-crossing with present base pairs are valid neighbored
insertion moves.
Valid deletion moves are all present base pairs with negative sign.
Valid shift moves are constructed by taking all paired positions as fix position of a shift move and iterating over
all positions of the sequence. If the letters of a position are compatible and if it the move is non-crossing with
existing base pairs, we have a valid shift move.
The method of generating shift moves can be accelerated by skipping neighbored base pairs.

If we need to construct all neighbors several times for subsequent moves, we can speed up the task by using
the move set of the previous structure. The previous move set has to be filtered, such that all moves that would
cross the next selected move are non-crossing. Next, the selected move has to be removed. Then one has to only
to generate all moves that were not possible before.
One move is the inverted selected move (if it was an insertion, simply make the indices negative).
The generation of all other new moves is different and depends on the selected move. It is easy for an insertion move,
because we have only to include all non-crossing shift moves, that are possible with the new base pair. For that we can
either iterate over the sequence or we can select all crossing shift moves in the filter procedure and convert them into
shifts.

The generation of new moves given a deletion is a little bit more complex, because we can create more moves. At first
we can insert the deleted pair as insertion move. Then we generate all insertions that would have crossed the deleted
base pair. Finally we construct all crossing shift moves.

If the given move is a shift, we can save much time by specifying the intervals for the generation of new moves.
The interval which was enclosed by the positive position of the shift move and the previous paired position is the
freed interval after applying the move. This freed interval includes all positions and base pairs that we need to
construct new insertions and shifts. All these new moves have one position in the freed interval and the other position
in the environment of the freed interval. The environment are all position which are outside the freed interval, but within
the same enclosing loop of the shift move. The environment for valid base pairs can be divided into one or more intervals,
depending on the shift move. The following examples describe a few scenarios to specify the intervals of the environment.

Example 1::

  increase the freed interval with a shift move
  AAAAGACAAGAAACAAAAGAGAAACAACAAACAAGAAACAAACAAAA
  ....(....(...)....(.(...)..)......(...)...).... // structure before the shift
  ....(....(...)....(.(...)......)..(...)...).... // structure after the shift
  ............................[__]............... // freed interval
  ..................[________]................... // interval that can pair with the freed interval

Example 2::

  switch the freed interval with a shift move
  AAAAGACAAGAAACAAAAGAGAAACAACAAACAAGAAACAAACAAAA
  ....(....(...)....(.(...)..)......(...)...).... // structure before the shift
  ....(.(..(...)....).(...).........(...)...).... // structure after the shift
  ...................[_______]................... // freed interval
  ....[_].....................[_____________].... // intervals that can pair with the freed interval

Example 3::

  decrease the freed interval with a shift move
  AAAAGACAAGAAACAAAAGAGAAACAACAAACAAGAAACAAACAAAA
  ....(....(...)....(.(...)......)..(...)...).... // structure before the shift
  ....(....(...)....(.(...)..)......(...)...).... // structure after the shift
  ............................[__]............... // freed interval
  ....[____________].............[__________].... // intervals that can pair with the freed interval

.. image:: /gfx/shift_move_intervals.svg

Given the intervals of the environment and the freed interval, the new shift moves can be constructed quickly.
One has to take all positions of pairs from the environment in order to create valid pairs with positions in the freed interval.
The same procedure can be applied for the other direction. This is taking all paired positions within the freed interval
in order to look for pairs with valid positions in the intervals of the environment.

.. doxygengroup:: neighbors
    :no-title:

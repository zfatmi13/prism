/*
 * Copyright (C) 2020  Zainab Fatmi and Franck van Breugel
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You can find a copy of the GNU General Public License at
 * <http://www.gnu.org/licenses/>.
 */

package explicit;

import java.util.HashMap;
import java.util.Map;

import prism.Evaluator;

/**
 * A non-empty sub-block. Each state of a sub-block has the same lifting. This
 * lifting of a state s maps each block b to the probability of the state s
 * transitioning to states in that block b. Only blocks that have been split in
 * the previous refinement step are considered.
 * 
 * @author Zainab Fatmi
 * @author Franck van Breugel
 */
public final class NonEmptySubBlock<Value> extends SubBlock<Value> {

	/**
	 * Evaluator for manipulating probability values in the distribution (of type {@code Value})
	 */
	private Evaluator<Value> eval;

	/**
	 * The lifting of this non-empty sub-block. It maps block IDs to probabilities.
	 */
	private final Map<Integer, Value> lifting;

	/**
	 * Initializes this non-empty sub-block.
	 */
	public NonEmptySubBlock(Evaluator<Value> eval) {
		super();
		this.lifting = new HashMap<Integer, Value>();
		this.eval = eval;
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public Value get(Integer block) {
		return this.lifting.get(block);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public void put(Integer block, Value probability) {
		this.lifting.put(block, probability);
	}

	/**
	 * Tests whether the given object is equal to this non-empty sub-block. Two
	 * non-empty sub-blocks are considered equal if their liftings maps blocks to
	 * probabilities that are very close.
	 * 
	 * @param object an object
	 * @return true if the given object is equal to this empty sub-block, false
	 *         otherwise
	 * @throws NullPointerException if object is null
	 * @throws ClassCastException   if object is not a sub-block
	 */
	@Override
	public boolean equals(Object object) {
		if (object instanceof EmptySubBlock) {
			return false;
		} else {
			@SuppressWarnings("unchecked")
			NonEmptySubBlock<Value> other = (NonEmptySubBlock<Value>) object;
			if (this.lifting.size() != other.lifting.size()) {
				return false;
			} else {
				for (Map.Entry<Integer, Value> entry : this.lifting.entrySet()) {
					Value probability = other.lifting.get(entry.getKey());
					if (!eval.equals(entry.getValue(), probability)) {
						return false;
					}
				}
			}
			return true;
		}
	}
}

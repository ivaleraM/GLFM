function hidden = remove_dims(hidden, idxs)

hidden.Z(:, idxs) = [];
hidden.B(:,idxs, :) = [];

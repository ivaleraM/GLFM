function hidden = sort_hidden(hidden)

    [xx, I] = sort(sum(hidden.Z),'descend');
    hidden.Z = hidden.Z(:,I);
    hidden.B = hidden.B(:,I,:);
    
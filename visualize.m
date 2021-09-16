function [] = visualize(state)
stateSize = size(state);
if stateSize(1) == 2
    qAvec = state(1,:);
    qBvec = state(1,:);
elseif stateSize(1) == 4
    qAvec = state(1,:);
    qBvec = state(2,:);
end

%params = getparams();
L = 1;
pAxvec = L*cos(qAvec);
pAyvec = L*sin(qAvec);

pBxvec = pAxvec + L*cos(qAvec+qBvec);
pByvec = pAyvec + L*sin(qAvec+qBvec);
fig = figure;
bodyA = plot([0,pAxvec(1)],[0,pAyvec(1)]);
hold on
bodyB = plot([pAxvec(1),pBxvec(2)],[pAyvec(1),pByvec(1)]);
axis([-2, 2,-2, 2]);
for iCnt = 1:stateSize(2)
    set(bodyA,'XData',[0,pAxvec(iCnt)])
    set(bodyA,'YData',[0,pAyvec(iCnt)])
    set(bodyB,'XData',[pAxvec(iCnt),pBxvec(iCnt)])
    set(bodyB,'YData',[pAyvec(iCnt),pByvec(iCnt)])
    drawnow
    frame = getframe(fig);
    im{iCnt} = frame2im(frame);
end


filename = 'testGif.gif';

for iCnt = 1:stateSize(2)
   [A,map] = rgb2ind(im{iCnt},256);
   if iCnt == 1
       imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0);
   else
       imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0);
   end
end